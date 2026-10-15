* NOTES :
    - main.cpp is now in drafts since the the code had been split into multiple files and the main loop is not yet complete.


* Done:
         - V3, M3 and quaternion implemented.
         - Quaternion rotation matrix, derivative and Euler-angle conversion implemented.
         - Rocket parameters now use const getter methods for derived geometry and CoGs.
         - Ogive shell, internal payload/mousse, flight computer, motor and fuel are
             represented as separate mass contributors in the CoG calculation.
         - ISA density model implemented through the 47-51 km layer using geopotential altitude.
         - Still-air drag law implemented.
         - RK4 couples position, velocity, orientation, angular velocity and fuel mass.
         - Constant fuel mass-flow rate is included in the RK4 stages.
         - Orientation has one source of truth: s.orientation.
         - A terminal simulation loop and basic state output exist.
         - Test R * R^T = I for the rotation matrix .
         - Separate math, rocket physics, integrator and main into different files.
         - Refresh s.drag after every integration step if it is kept as stored state/output.

 * Missing physics:
         - Wind model and relative air velocity.
         - Angle of attack from body-axis and relative-air-velocity directions.
         - Lift, normal force, aerodynamic side force and aerodynamic moments.
         - Mach/Reynolds-dependent aerodynamic coefficients and a Cd lookup/model.
         - Pitch/yaw damping and angle-of-attack-dependent center of pressure.
         - Engine model: thrust curve, mass-flow curve, burn duration, shutdown and thrust direction.
         - Remaining ISA layers above 51 km, if higher-altitude flight is required but not urgent.
         - Consistent ogive-shell inertia and internal payload/mousse inertia validation.
         - dI/dt effects during propellant burn.

 * Required corrections before trusting the simulation:
         - Initialize s.mass_fuel and s.dm from the rocket parameters.
         - Clamp fuel mass and burn rate when the remaining fuel is less than dm * dt.
         - Refresh s.drag after every integration step if it is kept as stored state/output.
         - Make inverse() return safely when the inertia matrix is singular.
         - Guard quaternion normalization against a zero or non-finite norm.
         - Add finite-value, positive-mass and valid-inertia state checks after every step.
         - Add ground-contact handling so position.z cannot remain below zero.

 * Validation:
         - Test rotation matrices with identity and known 90-degree rotations.
         - Test free fall with thrust and drag disabled.
         - Test constant-thrust vertical flight.
         - Test constant-burn fuel depletion against m(t) = m0 - dm * t.
         - Test conservation of energy/angular momentum where applicable.

 * Later:
         - Export trajectory data to CSV and plot it before adding a GUI/OpenGL layer.


* Last updated: 2026-10-16 *