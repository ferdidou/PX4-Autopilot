# Differential (jerk-level) control allocation — design notes

Method id: `CA_METHOD = 3` (`AllocationMethod::DIFFERENTIAL`).
Class: `ControlAllocationDifferential` (inherits `ControlAllocationPseudoInverse`).

This document is the mathematical reference for the dynamic allocator. It is built
in two layers: a linear ADA fallback and the nonlinear geometric NDA used by the
omnidirectional tiltrotor airframe.

## Why a *dynamic* allocator

PX4's stock allocators are **static / memoryless**: every cycle they solve
`a = a_trim + E^‡ w_d` from scratch (`E` = effectiveness matrix, `E^‡` its
pseudo-inverse, `w_d` the desired wrench, `a` the actuator vector).

The differential (jerk-level) allocation used for over-actuated tiltrotors
(Allenspach 2020, Cuniato 2024 "Augmented Differential Allocation", ADA) is
**dynamic**: it keeps the actuator state `a`, produces actuator *rates* `ȧ`, and
integrates them. The Augmented part (ADA) is what lets it run inside PX4 without
any acceleration feedback: it *synthesizes* the required wrench-rate (jerk)
directly from the wrench command PX4 already publishes:

    ẇ = k_j (w_d − w(a))

## Layer 1 — Linear ADA  (IMPLEMENTED)

Uses PX4's existing *linear* effectiveness `E ∈ R^(6×N)` as a constant Jacobian.

State:        a ∈ R^N        (= `_actuator_sp`, persists across cycles)
Wrench model: w(a) = E (a − a_trim)
Setpoint:     w_d = (control_sp − control_trim)
Pseudo-inv:   E^‡ = geninv(E) = `_mix` (N×6, normalized by the parent class)

Per-cycle update (`allocate()`):

    1.  w   = E (a − a_trim)                  // wrench from the maintained state
    2.  ẇ   = k_j (w_d − w)                    // ADA augmentation -> desired jerk
    3.  ȧ   = E^‡ ẇ                            // actuator rate
    4.  a  ← a + ȧ · Δt                        // INTEGRATION (Δt = measured dt)

### Closed-loop behaviour

Since w = E(a − a_trim):

    ẇ = E ȧ = k_j (E E^‡)(w_d − w)

For full row-rank E, E E^‡ = I_6, hence the decoupled linear ODE

    ẇ = k_j (w_d − w)   =>   w(t) = w_d + (w_0 − w_d) e^(−k_j t)

Exponential convergence, time constant τ = 1/k_j, no overshoot. `k_j` (param
`CA_ADA_KJ`) is the single tuning knob.

### Non-regression

At equilibrium (ẇ = 0): a → a_trim + E^‡ w_d, i.e. exactly the static
pseudo-inverse solution PX4 computes today. Layer 1 therefore reaches the same
steady state; it only adds a smooth, rate-limited transient and a persistent
internal state — the substrate Layer 2 needs.

### Discrete-time stability

    a_{k+1} = a_k + k_j Δt · E^‡ (w_d − E(a_k − a_trim))

Stable iff k_j Δt < 2 (eigenvalues of (I − k_j Δt E E^‡) inside the unit circle).
With Δt ≤ 0.02 s (CA loop is clamped to [0.2, 20] ms) and k_j ≈ 20,
k_j Δt ≤ 0.4 — comfortable margin.

## Layer 2 — Nonlinear NDA  (IMPLEMENTED, math core)

Replaces the constant linear Jacobian `J = E` by the configuration-dependent
analytic Jacobian of the lifted tiltrotor model, expressed directly in PX4's
*normalized* actuator coordinates (Normalized Differential Allocation, NDA — this
is what removes the manual [-1,1] / weight tuning).

State s = [s_alpha_0..N-1, s_thrust_0..N-1], where the tilt-servo output
s_alpha in [-1, 1] and the motor output s_thrust in [0, 1] are exactly the
values PX4 publishes on actuator_servos / actuator_motors. Per arm:

    alpha_i = c_alpha + h_alpha s_alpha_i   (c, h from CA_SV_TL*_MINA/MAXA)
    T_i     = thrust_max_i s_thrust_i       (thrust_max = CA_ROTOR*_CT)
    x_i = T_i cos(alpha_i),  y_i = T_i sin(alpha_i)

Per-arm thrust direction follows PX4's tiltedAxis:
axis(alpha) = cos(alpha) b_hat + sin(alpha) h_hat, with b_hat the nominal rotor
axis and h_hat = [cos(td), sin(td), 0] from CA_SV_TL*_TD. Constant wrench columns
(PX4 order: tau, then F; kappa = CA_ROTOR*_KM carries the spin sign):

    B^x_i = [ p_i x b_hat - kappa b_hat ; b_hat ]
    B^y_i = [ p_i x h_hat - kappa h_hat ; h_hat ]
    w = sum_i ( B^x_i x_i + B^y_i y_i )

(the -kappa drag/reaction term matches PX4's rotor convention
moment = p x thrust - km*thrust, with kappa = CA_ROTOR*_KM carrying the spin sign.)

Analytic Jacobian (chain rule through the normalized->physical map):

    dw/ds_alpha_i  = h_alpha (-y_i B^x_i + x_i B^y_i)
    dw/ds_thrust_i = thrust_max (cos(alpha_i) B^x_i + sin(alpha_i) B^y_i)

Allocation loop (allocate(), geometric branch):

    1-2. w, J = computeWrenchAndJacobian(s)
    3.   wdot = k_j (w_d - w)
    4.   s_dot = J^+ wdot                         (geninv, right pseudo-inverse)
    5.   s_dot += (I - J^+ J) s_dot*              (nullspace: tilts->level, motors->hover)
    6.   s += s_dot dt ; clip s_alpha in[-1,1], s_thrust in[0,1] ; write actuator_sp

Verified by ControlAllocationDifferentialTest.cpp: (a) the analytic Jacobian
matches a central finite difference to 1e-2, (b) the integrated state drives the
wrench to a feasible setpoint (residual < 1e-2). The differential form is
intrinsically robust to the +-90 deg tilt singularities that break the static
geometric inverse.

### Quadratic-motor output

The geometric model is linear in the state, so `s_thrust` is the desired thrust
fraction. The simulator/ESC produce thrust ~ cmd^2, so the motor command is emitted
as `sqrt(s_thrust)`: `thrust_max (sqrt s)^2 = thrust_max s` recovers the commanded
thrust. The sqrt is applied only to the output (never differentiated), so the loop —
which lives in s-space — is unchanged and has no Jacobian singularity at s = 0.

### Flight wiring

`ActuatorEffectivenessOmniTilt` (CA_AIRFRAME=16) fills ArmGeometry[] from
CA_ROTOR* / CA_SV_TL* and `ControlAllocator` calls `setOmniGeometry()`. The wrench
columns use the same normalized scaling as PX4's thrust/torque setpoints
(thrust_max = CA_ROTOR*_CT). Airframe + SITL model: 6 rotors + 6 tilt servos
(`10020_gazebo-classic_hexa_tilting`).

Later ADA-paper refinement: PDA (propeller power dynamics).
