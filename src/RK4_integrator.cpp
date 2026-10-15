#include "RK4_integrator.hpp"

void rk4(rocket_dynamics& s, float dt){
        
        // initialization of the RK4 stages
        float fuel_mass_init = s.mass_fuel;
        const float fuel_dm = (s.dm > 0.0f) ? -s.dm : s.dm; // If dm is stored as a positive burn magnitude, make it negative so fuel decreases.
        V3 w0 = s.angular_velocity;
        quaternion q0 = s.orientation;
        V3 init_stage_p = s.position;
        V3 init_stage_vel = s.velocity;
        
        
        // Stage 1
        float fuel_mass_stage = fuel_mass_init;
        quaternion stage_q = q0;
        V3 stage_w = w0;
        V3 stage = init_stage_p;
        V3 stage_vel = init_stage_vel;
        V3 stage_drag = s.drag_force(stage_vel, stage);

        float km1 = fuel_dm;
        V3 Kv1 = s.accel(stage_q, stage_drag, fuel_mass_stage);
        V3 Kp1 = stage_vel;
        quaternion kq1 = stage_q.deriv(stage_w);
        V3 kw1 = s.angular_accel(stage_w, stage_q, stage_drag, fuel_mass_stage);

        
        //stage 2
        fuel_mass_stage = fuel_mass_init + km1 * dt / 2.0f;
        stage_q = q0 + kq1 * dt / 2.0f;
        stage_q.normalize(); // Normalize the quaternion to prevent drift
        stage_w = w0 + kw1 * dt / 2.0f;
        stage = init_stage_p + Kp1 * dt / 2.0f;
        stage_vel = init_stage_vel + Kv1 * dt / 2.0f;
        stage_drag = s.drag_force(stage_vel, stage);

        float km2 = fuel_dm;
        V3 Kv2 = s.accel(stage_q, stage_drag, fuel_mass_stage);
        V3 Kp2 = stage_vel;
        quaternion kq2 = stage_q.deriv(stage_w);
        V3 kw2 = s.angular_accel(stage_w, stage_q, stage_drag, fuel_mass_stage);
        

        //stage 3
        fuel_mass_stage = fuel_mass_init + 0.5f * km2 * dt;
        stage_q = q0 + kq2 * dt / 2.0f;
        stage_q.normalize(); // Normalize the quaternion to prevent drift
        stage_w = w0 + kw2 * dt / 2.0f;
        stage = init_stage_p + Kp2 * dt / 2.0f;
        stage_vel = init_stage_vel + Kv2 * dt / 2.0f;
        stage_drag = s.drag_force(stage_vel, stage);
        
        float km3 = fuel_dm;
        V3 Kv3 = s.accel(stage_q, stage_drag, fuel_mass_stage);
        V3 Kp3 = stage_vel;
        quaternion kq3 = stage_q.deriv(stage_w);
        V3 kw3 = s.angular_accel(stage_w, stage_q, stage_drag, fuel_mass_stage);
        

        //stage 4
        fuel_mass_stage = fuel_mass_init + km3 * dt;
        stage_q = q0 + kq3 * dt;
        stage_q.normalize(); // Normalize the quaternion to prevent drift
        stage_w = w0 + kw3 * dt;
        stage = init_stage_p + Kp3 * dt;
        stage_vel = init_stage_vel + Kv3 * dt;
        stage_drag = s.drag_force(stage_vel, stage);

        float km4 = fuel_dm;
        V3 Kv4 = s.accel(stage_q, stage_drag, fuel_mass_stage);
        V3 Kp4 = stage_vel;
        quaternion kq4 = stage_q.deriv(stage_w);
        V3 kw4 = s.angular_accel(stage_w, stage_q, stage_drag, fuel_mass_stage);
        

        // Final update
        s.mass_fuel = fuel_mass_init + (km1 + 2.0f * km2 + 2.0f * km3 + km4) * (dt / 6.0f);
        s.orientation = s.orientation+ (kq1 + kq2 * 2.0f + kq3 * 2.0f + kq4) * (dt / 6.0f);
        s.angular_velocity = s.angular_velocity + (kw1 + kw2 * 2.0f + kw3 * 2.0f + kw4) * (dt / 6.0f); 
        s.velocity = s.velocity + (Kv1 + Kv2*2.0 + Kv3*2.0 + Kv4) * (dt / 6.0f);
        s.position = s.position + (Kp1 + Kp2*2 + Kp3*2 + Kp4) * (dt/6.0f);
        s.orientation.normalize();
    }