# Progress Log

## before 4/29

1. set up a linux development environment on my laptop and connect it with my macmini.
2. review basic cpp syntax and learned the basic cmake workflow for building and running a cpp project.
3. started a simple 2D heat diffusion solver project and implemented several key functions:

- idx : maps a 2D grid location (i,j) to a 1D memory index. This is used to store a 2D field in a 1D std::vector<double>.
- intial_t_field : initializes the temperature field. the field is intialized to 0.0 everywhere, and center square reigon is set to 40.0.
- step_heat : performs one explicit heat diffusion step using a 5-point stencil. This function updates the temperature field from the old time step to the new time step. One important observation is that the time step parameter `alpha` must stay below the stability limit; otherwise, the solution becomes unstable.

4. Added CSV output and Python visualization scripts to inspect the simulation results.

5. Observed numerical instability when `alpha` is too large. In particular, using `alpha = 0.3` produced a checkerboard-like pattern, which indicates that the explicit scheme is unstable.

