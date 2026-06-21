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
 * @file ControlAllocationDifferential.hpp
 *
 * Dynamic (jerk-level) control allocation, a.k.a. Augmented Differential
 * Allocation (ADA). Unlike the static pseudo-inverse, this allocator keeps the
 * actuator state, produces actuator *rates* and integrates them. See
 * DifferentialAllocation.md for the full derivation.
 *
 * Two modes:
 *  - Linear (Layer 1): reuses PX4's linear effectiveness E as a constant
 *    Jacobian. Active when no tiltrotor geometry has been set.
 *  - Geometric / NDA (Layer 2): nonlinear lifted tiltrotor model expressed
 *    directly in PX4's *normalized* actuator coordinates (Normalized Differential
 *    Allocation). The state is s = [s_alpha_0.., s_thrust_0..] with the tilt-servo
 *    output s_alpha in [-1, 1] (maps linearly to [alpha_min, alpha_max]) and the
 *    motor output s_thrust in [0, 1] (normalized thrust). Active once
 *    setOmniGeometry() has been called.
 *
 * Lifted coordinates x_i = T_i cos(alpha_i), y_i = T_i sin(alpha_i) with
 * T_i = thrust_max_i * s_thrust_i and alpha_i = c_alpha_i + h_alpha_i * s_alpha_i.
 */

#pragma once

#include "ControlAllocationPseudoInverse.hpp"

#include <px4_platform_common/module_params.h>

class ControlAllocationDifferential : public ControlAllocationPseudoInverse, public ModuleParams
{
public:
	ControlAllocationDifferential() : ModuleParams(nullptr) {}
	virtual ~ControlAllocationDifferential() = default;

	void allocate() override;

	void updateParameters() override { updateParams(); }

	bool isDifferential() const override { return true; }

	// Report the actually produced wrench (geometric model) so the unallocated-control
	// / anti-windup feedback to the rate controller is correct (the base version uses
	// the linear effectiveness matrix, which does not match the geometric allocation).
	matrix::Vector<float, NUM_AXES> getAllocatedControl() const override { return _allocated_wrench; }

	static constexpr int MAX_ARMS = 6;
	static constexpr int NUM_LIFTED = 2 * MAX_ARMS; ///< dimension of the [s_alpha; s_thrust] state

	/**
	 * Per-arm geometry of a tiltable rotor, in PX4-normalized actuator coordinates.
	 * The thrust direction sweeps as axis(alpha) = cos(alpha) b_hat + sin(alpha) h_hat
	 * (PX4 tiltedAxis convention).
	 */
	struct ArmGeometry {
		matrix::Vector3f b_hat;    ///< nominal (untilted) thrust axis, unit
		matrix::Vector3f h_hat;    ///< horizontal sweep direction (axis at +90 deg tilt), unit
		matrix::Vector3f position; ///< rotor position [m], body frame
		float thrust_max{1.f};     ///< thrust at s_thrust = 1 (PX4 thrust_coef)
		float kappa{0.f};          ///< signed drag-to-thrust ratio Cm/Ct (sign encodes spin)
		float alpha_min{0.f};      ///< tilt angle [rad] at s_alpha = -1
		float alpha_max{0.f};      ///< tilt angle [rad] at s_alpha = +1
		int alpha_index{-1};       ///< actuator-vector index of the tilt servo (output s_alpha)
		int thrust_index{-1};      ///< actuator-vector index of the motor (output s_thrust)
	};

	/**
	 * Provide the per-arm lifted geometry. Switches allocate() into the nonlinear
	 * NDA mode and initializes the state to level hover (s_alpha = 0, s_thrust = hover).
	 */
	void setOmniGeometry(const ArmGeometry arms[], int num_arms, float hover_thrust);

	bool hasGeometry() const { return _num_arms > 0; }

	/**
	 * Nonlinear tiltrotor model in normalized coordinates: wrench and analytic
	 * Jacobian at state (s_alpha, s_thrust). Exposed for unit testing / wiring.
	 *
	 * @param s_alpha  normalized tilt-servo outputs [-1, 1], length num_arms
	 * @param s_thrust normalized motor outputs [0, 1], length num_arms
	 * @param wrench   output wrench [tau_x tau_y tau_z F_x F_y F_z]
	 * @param jacobian output J = dw/ds, columns [s_alpha_0..N-1, s_thrust_0..N-1]
	 */
	void computeWrenchAndJacobian(const float s_alpha[], const float s_thrust[],
				      matrix::Vector<float, NUM_AXES> &wrench,
				      matrix::Matrix<float, NUM_AXES, NUM_LIFTED> &jacobian) const;

	void getState(float s_alpha_out[], float s_thrust_out[]) const;

private:
	void allocateLinear();    ///< Layer 1: linear ADA using PX4's effectiveness matrix
	void allocateGeometric(); ///< Layer 2: nonlinear NDA
	void buildLiftedColumns();

	DEFINE_PARAMETERS(
		(ParamFloat<px4::params::CA_ADA_KJ>) _param_ca_ada_kj   ///< wrench convergence gain k_j
	)

	// Geometry (constant lifted columns) and normalized state
	matrix::Vector<float, NUM_AXES> _Bx[MAX_ARMS]; ///< wrench column for lifted x_i (b_hat direction)
	matrix::Vector<float, NUM_AXES> _By[MAX_ARMS]; ///< wrench column for lifted y_i (h_hat direction)
	matrix::Vector<float, NUM_AXES> _wrench_scale; ///< per-axis scale mapping physical wrench -> PX4 normalized
	ArmGeometry _arms[MAX_ARMS];
	int _num_arms{0};
	float _s_alpha[MAX_ARMS] {};  ///< tilt-servo state [-1, 1]
	float _s_thrust[MAX_ARMS] {}; ///< motor state [0, 1]
	float _hover_thrust{0.f};
	matrix::Vector<float, NUM_AXES> _allocated_wrench; ///< wrench actually produced (for unallocated/anti-windup)

	// Secondary (nullspace) task gains: gently regularize the redundant DOF.
	static constexpr float NULLSPACE_ALPHA_GAIN = 1.0f; ///< pull tilt servos toward 0 (level)
	static constexpr float NULLSPACE_THRUST_GAIN = 1.0f; ///< pull motors toward hover
};
