# LQR-Trajectory-Generator
This repository contains Matlab scripts for the aim of generating 3D trajectory in the Special Euclidean SE(3) for quadrotors while minimizing the crackle (the 5th derivative of the position) using the Linear Quadratic Regulator **LQR**.

This repository tries to reproduce the simulations of the paper entitled **"Trajectory Generation of on SE(3) for an Underactuated Vehicle with Pointing Direction Constraints"** in https://ieeexplore.ieee.org/document/8815238 . Two methods are used: in `main.m`, the **Euler** method is used for the different integrations, however in `main_2.m`, the **ODE45** Matlab solver is used.

# How To Use The Code
1- Depending on the quadrotor you may run the tests on, you must change the variables of the mass and the inertial matrix in the code to the ones of your quad

2- Run `main.m` or `main_2.m`
