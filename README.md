# Numerical-LRE-Solution

## Description

Numerical solutions to the backreaction equation and linear response equation for the semiclassical approximation to spin-1/2 quantum electrodynamics in 1+1 dimensions.



## File structure

The `src` folder contains the bulk of the implementation.

​	`simp13` - Implementation of the Simpson's 1/3 method.

​	`SBE` - System of backreaction equations for A,E with a spin-1/2 field in 1+1 dimensions

​	`Solve_LRE` - Primary section for solving ODE for `SBE` and solving linear response equation.

The `runs` folder contains the entry point for the code.

​	`Start_Run.m` - Entry point and parameter list.



## Getting Started

To start a run, adjust parameters within and run from `runs/Start_Run.m`.

A headless run can be started via

```
matlab -nodesktop -nosplash -r "run('Start_Run.m');" 
```

