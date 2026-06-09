***
# 6DOF ROCKET SIMULATOR PLAN :

### This project will evolve slowly as I progress in acquiring the necessary knowledge to model a rocket and satellite physics and their control systems with hardware in loop 

### While right now I'm still dealing with the physics part but the evolution of this project is thought, as I'm planing to have at least 5 main versions which would be released with corrections in between :
--- 
---
*** 
# [ ] V.0 : Raw physics only in terminal :
- This version will only feature physics part such as the translational and inertial equations which would describe the rockets trajectory while reacting to only 3 main forces : thrust, drag and its weight.

- While the results would be shown only in the terminal, the goal would be to model the physics of a rocket correctly before things get complicated .


### Notes : 
- It's highly likely the model of the motor is going to be pretty accurate and not just generic (A class L motor or K might be modeled) .

- I'll add links to the actual rocket's model in a cad and will also add some pictures in the readme file

# [ ] V.1 : Open GL implementation :
- Having the physics working is great but seeing how our model actually reacts in real life is best , thus I'll dive into the 3D visualisation of the object and the environment so we become more familiar with the physics since the next parts of the project requires quite some intuition and an actual moving object to understand what's happening .

### Note: 
- I preferred to keep that part separate since OpenGL is not that simple and it actually requires some knowledge about my hardware and how it works .
model. 
# [ ] V.2 : Control implementation :

- Here things will get kind of coarse since it requires a study of control systems of rockets, unlike a drone's control system.

- Heavy control theory concepts are going to be tested and implemented from a simple linear control system (PID) to a robust one such as a multivariable non linear control system .

- I'll only implement the Guidance navigation and control part of the rocket before reaching the space (attitude control only) .

# [ ] V.3 : Orbital Mechanics and orbital control implementation :

- This part will be the heaviest one since we'll actually do two parts in one :
    
    - ### Part A : Orbital mechanics implementation : 
        - Applying Kepler's laws and the 2 body problem to get an accurate estimation of the trajectory of the rocket is going to be the plan here 
    
    - ### Part B : Orbital control implementation :
        - I'll start with this part only when the part A is complete and correct or this part would be an absolute mess .

        - This part I still do not know exactly what I'll face since I only have a generic idea about those systems but that's what makes the project the most interesting .

- # Notes :
    - You can notice that this part is the most exciting one since we'll dive into actual rocket science and astronautics engineering .

# V.4 : THE MISSION :

- Our rocket cannot fly without any purpose so here we'll simulate a mission of a satellite , the goal is to estimate its trajectory and know how much it could stay up there 

- This part is going to require to model a whole space engine and simulate its physics too , so I might simulate the satellite alone first then I'll proceed to add it on with the rocket.

---

# This project is Super ambitious, and one of its ambitions is making the project accessible so other people could also simulate trajectories without having to deal with the headache of programming and the physics if those aren't their thing.



