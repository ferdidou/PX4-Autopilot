/****************************************************************************
 *
 *   Copyright (c) 2024 PX4 Development Team. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 * 3. Neither the name PX4 nor the names of its contributors may be
 *    used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 * COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
 * OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 * AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 * ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 ****************************************************************************/

/**
 * @file ControlAllocationDifferential.cpp
 *
 * Dynamic (jerk-level) control allocation.
 *  - Layer 1: linear ADA (uses PX4's effectiveness matrix as a constant Jacobian)
 *  - Layer 2: nonlinear NDA (lifted tiltrotor model in normalized actuator coordinates)
 */

#include "ControlAllocationDifferential.hpp"

using matrix::Vector3f;

void
ControlAllocationDifferential::allocate()
{
	if (hasGeometry()) {
		allocateGeometric();

	} else {
		allocateLinear();
	}
}

void
ControlAllocationDifferential::allocateLinear()
{
	// Layer 1: linear ADA. w(a) = E (a - a_trim), wdot = k_j (w_d - w), a += E^+ wdot dt.
	// In steady state this reaches the same point as the static pseudo-inverse.
	updatePseudoInverse();

	_prev_actuator_sp = _actuator_sp;

	const float k_j = _param_ca_ada_kj.get();
	const matrix::Vector<float, NUM_AXES> wrench = _effectiveness * (_actuator_sp - _actuator_trim);
	const matrix::Vector<float, NUM_AXES> wrench_dot = k_j * ((_control_sp - _control_trim) - wrench);
	const ActuatorVector actuator_sp_dot = _mix * wrench_dot;
	_actuator_sp = _actuator_sp + actuator_sp_dot * _dt;

	_allocated_wrench = (_effectiveness * (_actuator_sp - _actuator_trim)).emult(_control_allocation_scale);
}

void
ControlAllocationDifferential::buildLiftedColumns()
{
	// Constant per-arm wrench columns of the lifted model (PX4 axis order: tau, then F).
	// Matching PX4's rotor moment convention (ActuatorEffectivenessRotors):
	//   moment = p x thrust - km * thrust   (note the MINUS on the reaction/drag term)
	//   B^x_i = [ p_i x b_hat - kappa*b_hat ; b_hat ]
	//   B^y_i = [ p_i x h_hat - kappa*h_hat ; h_hat ]
	// (kappa already carries the rotor-spin sign.)
	for (int i = 0; i < _num_arms; ++i) {
		const Vector3f b = _arms[i].b_hat;
		const Vector3f h = _arms[i].h_hat;
		const Vector3f p = _arms[i].position;
		const float k = _arms[i].kappa;

		const Vector3f tau_b = p.cross(b) - b * k;
		const Vector3f tau_h = p.cross(h) - h * k;

		for (int j = 0; j < 3; ++j) {
			_Bx[i](j) = tau_b(j);
			_Bx[i](j + 3) = b(j);
			_By[i](j) = tau_h(j);
			_By[i](j + 3) = h(j);
		}
	}

	// Per-axis normalization so the geometric wrench lives in the same normalized
	// range as PX4's thrust/torque setpoints (pitfall: PX4 setpoints are NOT in
	// N / N.m). Scale of axis r = sum over arms of the max physical contribution
	// (thrust_max * |column entry|). A fully-deflected actuator set then maps an
	// axis to ~[-1, 1], matching the controller output.
	_wrench_scale.setAll(1.f);

	for (int r = 0; r < NUM_AXES; ++r) {
		float scale = 0.f;

		for (int i = 0; i < _num_arms; ++i) {
			scale += _arms[i].thrust_max * (fabsf(_Bx[i](r)) + fabsf(_By[i](r)));
		}

		if (scale > 1e-6f) {
			_wrench_scale(r) = scale;
		}
	}
}

void
ControlAllocationDifferential::computeWrenchAndJacobian(const float s_alpha[], const float s_thrust[],
		matrix::Vector<float, NUM_AXES> &wrench,
		matrix::Matrix<float, NUM_AXES, NUM_LIFTED> &jacobian) const
{
	wrench.setZero();
	jacobian.setZero();

	for (int i = 0; i < _num_arms; ++i) {
		// Normalized -> physical: alpha = c + h s_alpha, T = thrust_max s_thrust
		const float h_alpha = 0.5f * (_arms[i].alpha_max - _arms[i].alpha_min);
		const float c_alpha = 0.5f * (_arms[i].alpha_max + _arms[i].alpha_min);
		const float alpha = c_alpha + h_alpha * s_alpha[i];
		const float thrust = _arms[i].thrust_max * s_thrust[i];

		const float ca = cosf(alpha);
		const float sa = sinf(alpha);
		const float x = thrust * ca;   // lifted x_i = T cos(alpha)
		const float y = thrust * sa;   // lifted y_i = T sin(alpha)

		// Wrench: w += B^x_i x_i + B^y_i y_i
		wrench += _Bx[i] * x + _By[i] * y;

		// Jacobian columns (chain rule through the normalized->physical map):
		//   dw/ds_alpha_i  = h_alpha (-y_i B^x_i + x_i B^y_i)
		//   dw/ds_thrust_i = thrust_max (cos B^x_i + sin B^y_i)
		const matrix::Vector<float, NUM_AXES> col_alpha = (_Bx[i] * (-y) + _By[i] * x) * h_alpha;
		const matrix::Vector<float, NUM_AXES> col_thrust = (_Bx[i] * ca + _By[i] * sa) * _arms[i].thrust_max;

		for (int r = 0; r < NUM_AXES; ++r) {
			jacobian(r, i) = col_alpha(r);              // columns 0..N-1: tilt servos
			jacobian(r, MAX_ARMS + i) = col_thrust(r);  // columns MAX_ARMS..: motors
		}
	}

	// Express wrench and Jacobian in PX4-normalized axes (divide each output row
	// by its scale). w_d arrives normalized, so this keeps the ADA loop consistent.
	for (int r = 0; r < NUM_AXES; ++r) {
		const float inv = 1.f / _wrench_scale(r);
		wrench(r) *= inv;

		for (int c = 0; c < NUM_LIFTED; ++c) {
			jacobian(r, c) *= inv;
		}
	}
}

void
ControlAllocationDifferential::allocateGeometric()
{
	const float k_j = _param_ca_ada_kj.get();

	// 1-2. Wrench and analytic Jacobian at the current maintained (normalized) state
	matrix::Vector<float, NUM_AXES> wrench;
	matrix::Matrix<float, NUM_AXES, NUM_LIFTED> jacobian;
	computeWrenchAndJacobian(_s_alpha, _s_thrust, wrench, jacobian);

	// 3. ADA augmentation: synthesize the desired wrench-rate from the wrench error.
	const matrix::Vector<float, NUM_AXES> wrench_dot = k_j * (_control_sp - wrench);

	// 4. Map wrench-rate to state-rate through the (right) pseudo-inverse of the Jacobian.
	matrix::Matrix<float, NUM_LIFTED, NUM_AXES> jacobian_pinv;
	matrix::geninv(jacobian, jacobian_pinv);
	matrix::Vector<float, NUM_LIFTED> s_dot = jacobian_pinv * wrench_dot;

	// 5. Nullspace regularization: gently pull the redundant DOF toward a nominal
	//    posture (tilts -> level, motors -> hover) without disturbing the wrench.
	//    s_dot += (I - J^+ J) s_dot*
	matrix::Vector<float, NUM_LIFTED> s_dot_secondary;

	for (int i = 0; i < _num_arms; ++i) {
		s_dot_secondary(i) = -NULLSPACE_ALPHA_GAIN * _s_alpha[i];
		s_dot_secondary(MAX_ARMS + i) = NULLSPACE_THRUST_GAIN * (_hover_thrust - _s_thrust[i]);
	}

	const matrix::SquareMatrix<float, NUM_LIFTED> nullspace =
		matrix::eye<float, NUM_LIFTED>() - jacobian_pinv * jacobian;
	s_dot += nullspace * s_dot_secondary;

	// 6. Integrate the normalized state, clip to actuator ranges, write to outputs.
	_prev_actuator_sp = _actuator_sp;

	for (int i = 0; i < _num_arms; ++i) {
		_s_alpha[i] += s_dot(i) * _dt;
		_s_thrust[i] += s_dot(MAX_ARMS + i) * _dt;

		// s_alpha in [-1, 1] (tilt servo), s_thrust in [0, 1] (motor)
		_s_alpha[i] = (_s_alpha[i] < -1.f) ? -1.f : (_s_alpha[i] > 1.f) ? 1.f : _s_alpha[i];
		_s_thrust[i] = (_s_thrust[i] < 0.f) ? 0.f : (_s_thrust[i] > 1.f) ? 1.f : _s_thrust[i];

		if (_arms[i].alpha_index >= 0) {
			_actuator_sp(_arms[i].alpha_index) = _s_alpha[i];
		}

		if (_arms[i].thrust_index >= 0) {
			// Quadratic-motor inversion. The geometric model treats thrust as
			// linear in the state (thrust = thrust_max * s_thrust), so s_thrust IS
			// the desired thrust fraction. The simulator/ESC produce thrust ~ cmd^2,
			// so emit cmd = sqrt(s_thrust): then thrust_max*(sqrt(s))^2 = thrust_max*s
			// is exactly the commanded thrust. This keeps the whole differential loop
			// (which lives in s-space) unchanged -- only the final emission is mapped,
			// so there is no Jacobian singularity at s = 0 (sqrt is applied to the
			// output, never differentiated). Matches the PRISMA allocator (cmd = sqrt(T)).
			_actuator_sp(_arms[i].thrust_index) = sqrtf(_s_thrust[i]);
		}
	}

	// Wrench actually produced by the new (clipped) state, reported via
	// getAllocatedControl() so the rate controller's anti-windup is correct.
	matrix::Matrix<float, NUM_AXES, NUM_LIFTED> jacobian_unused;
	computeWrenchAndJacobian(_s_alpha, _s_thrust, _allocated_wrench, jacobian_unused);
}

void
ControlAllocationDifferential::setOmniGeometry(const ArmGeometry arms[], int num_arms, float hover_thrust)
{
	_num_arms = (num_arms < MAX_ARMS) ? num_arms : MAX_ARMS;
	_hover_thrust = hover_thrust;

	for (int i = 0; i < _num_arms; ++i) {
		_arms[i] = arms[i];
		_arms[i].b_hat.normalize();
		_arms[i].h_hat.normalize();
		_s_alpha[i] = 0.f;
		_s_thrust[i] = hover_thrust;

		if (_arms[i].alpha_index >= 0) {
			_actuator_sp(_arms[i].alpha_index) = 0.f;
		}

		if (_arms[i].thrust_index >= 0) {
			// sqrt: same quadratic-motor inversion as the emitted command below.
			_actuator_sp(_arms[i].thrust_index) = sqrtf(hover_thrust);
		}
	}

	buildLiftedColumns();
}

void
ControlAllocationDifferential::getState(float s_alpha_out[], float s_thrust_out[]) const
{
	for (int i = 0; i < _num_arms; ++i) {
		s_alpha_out[i] = _s_alpha[i];
		s_thrust_out[i] = _s_thrust[i];
	}
}
