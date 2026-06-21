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
 * Test the nonlinear NDA (Layer 2) part of the differential allocator:
 *  - the analytic Jacobian matches a central finite difference of the wrench,
 *  - the integrated (normalized) state drives the wrench to a feasible setpoint.
 */

#include <gtest/gtest.h>
#include <math.h>
#include <ControlAllocationDifferential.hpp>

using namespace matrix;

namespace
{

// A small subclass to set a known geometry and reach the test helpers.
class TestableDifferential : public ControlAllocationDifferential
{
public:
	// Build a 6-arm omnidirectional hexarotor: arms every 60 deg, tilt sweeping in the
	// tangential-vertical plane, alternating rotor spin direction, +-30 deg tilt range.
	void buildHexOmni(float hover_thrust)
	{
		ArmGeometry arms[6];

		for (int i = 0; i < 6; ++i) {
			const float phi = (float)i * (float)M_PI / 3.f; // 60 deg spacing
			const float r = 0.3f;
			arms[i].position = Vector3f(r * cosf(phi), r * sinf(phi), 0.f);
			arms[i].b_hat = Vector3f(0.f, 0.f, -1.f);              // nominal thrust up
			arms[i].h_hat = Vector3f(-sinf(phi), cosf(phi), 0.f);  // tangential sweep
			arms[i].thrust_max = 10.f;
			arms[i].kappa = (i % 2 == 0) ? 0.05f : -0.05f;         // signed drag/thrust
			arms[i].alpha_min = -0.5236f;                          // -30 deg
			arms[i].alpha_max = 0.5236f;                           // +30 deg
			arms[i].alpha_index = i;
			arms[i].thrust_index = 6 + i;
		}

		setOmniGeometry(arms, 6, hover_thrust);
	}
};

} // namespace

TEST(ControlAllocationDifferentialTest, JacobianMatchesFiniteDifference)
{
	TestableDifferential method;
	method.buildHexOmni(0.5f);

	// Arbitrary non-trivial normalized state, within bounds.
	float s_alpha[6];
	float s_thrust[6];

	for (int i = 0; i < 6; ++i) {
		s_alpha[i] = 0.1f * (float)i - 0.3f;  // in [-0.3, 0.2]
		s_thrust[i] = 0.4f + 0.05f * (float)i; // in [0.4, 0.65]
	}

	Vector<float, 6> w0;
	Matrix<float, 6, 12> J_analytic;
	method.computeWrenchAndJacobian(s_alpha, s_thrust, w0, J_analytic);

	const float eps = 1e-4f;

	// Columns 0..5: d/ds_alpha_i, columns 6..11: d/ds_thrust_i
	for (int i = 0; i < 6; ++i) {
		for (int which = 0; which < 2; ++which) {
			float sa[6];
			float st[6];

			for (int k = 0; k < 6; ++k) { sa[k] = s_alpha[k]; st[k] = s_thrust[k]; }

			float *var = (which == 0) ? &sa[i] : &st[i];
			const int col = (which == 0) ? i : 6 + i;

			Vector<float, 6> w_plus;
			Vector<float, 6> w_minus;
			Matrix<float, 6, 12> dummy;

			const float saved = *var;
			*var = saved + eps;
			method.computeWrenchAndJacobian(sa, st, w_plus, dummy);
			*var = saved - eps;
			method.computeWrenchAndJacobian(sa, st, w_minus, dummy);
			*var = saved;

			const Vector<float, 6> fd = (w_plus - w_minus) / (2.f * eps);

			for (int r = 0; r < 6; ++r) {
				EXPECT_NEAR(J_analytic(r, col), fd(r), 1e-2f)
						<< "mismatch at row " << r << " col " << col;
			}
		}
	}
}

TEST(ControlAllocationDifferentialTest, WrenchConvergesToFeasibleSetpoint)
{
	TestableDifferential method;
	method.buildHexOmni(0.5f);
	method.updateParameters(); // load CA_ADA_KJ default

	// Build a feasible wrench setpoint = the wrench at a known target state (within bounds).
	float sa_t[6];
	float st_t[6];

	for (int i = 0; i < 6; ++i) {
		sa_t[i] = 0.2f * sinf((float)i);
		st_t[i] = 0.5f + 0.03f * (float)i;
	}

	Vector<float, 6> w_target;
	Matrix<float, 6, 12> dummy;
	method.computeWrenchAndJacobian(sa_t, st_t, w_target, dummy);

	method.setControlSetpoint(w_target);
	method.setDt(0.004f); // 250 Hz

	for (int step = 0; step < 2000; ++step) {
		method.allocate();
	}

	float sa_f[6];
	float st_f[6];
	method.getState(sa_f, st_f);

	Vector<float, 6> w_final;
	method.computeWrenchAndJacobian(sa_f, st_f, w_final, dummy);

	const Vector<float, 6> err = w_final - w_target;
	EXPECT_LT(err.norm(), 1e-2f) << "wrench did not converge, residual = " << (double)err.norm();
}
