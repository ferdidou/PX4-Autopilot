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
 * @file ActuatorEffectivenessOmniTilt.hpp
 *
 * Effectiveness for an omnidirectional multirotor where every rotor has its own
 * tilt servo (one rotor per tilt). It reuses ActuatorEffectivenessRotors and
 * ActuatorEffectivenessTilts to parse the standard CA_ROTOR* / CA_SV_TL* params,
 * and is meant to be driven by the dynamic differential allocator
 * (CA_METHOD = DIFFERENTIAL). The linear effectiveness matrix it builds is only
 * used as a fallback; the differential allocator reads the geometry directly via
 * the accessors below.
 */

#pragma once

#include "ActuatorEffectiveness.hpp"
#include "ActuatorEffectivenessRotors.hpp"
#include "ActuatorEffectivenessTilts.hpp"

#include <px4_platform_common/module_params.h>

class ActuatorEffectivenessOmniTilt : public ModuleParams, public ActuatorEffectiveness
{
public:
	ActuatorEffectivenessOmniTilt(ModuleParams *parent);
	virtual ~ActuatorEffectivenessOmniTilt() = default;

	bool getEffectivenessMatrix(Configuration &configuration, EffectivenessUpdateReason external_update) override;

	int numMatrices() const override { return 1; }

	void getDesiredAllocationMethod(AllocationMethod allocation_method_out[MAX_NUM_MATRICES]) const override
	{
		allocation_method_out[0] = AllocationMethod::DIFFERENTIAL;
	}

	void getNormalizeRPY(bool normalize[MAX_NUM_MATRICES]) const override
	{
		// The differential allocator does not use the normalized linear mix.
		normalize[0] = false;
	}

	const char *name() const override { return "Omni Tilt"; }

	// Accessors used by ControlAllocator to build the differential-allocator geometry.
	const ActuatorEffectivenessRotors::Geometry &rotorGeometry() const { return _rotors.geometry(); }
	const ActuatorEffectivenessTilts &tilts() const { return _tilts; }

protected:
	ActuatorEffectivenessRotors _rotors;
	ActuatorEffectivenessTilts _tilts;
	int _first_tilt_idx{0};
};
