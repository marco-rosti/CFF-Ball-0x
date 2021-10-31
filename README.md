# CFF-Ball-0x
CFF-Ball-0x is a code that simulates dense suspensions of rigid, spherical particles in a shear flow. The particles interact through hydrodynamics, contacts and conservative repulsive and attractive forces.

## About the code
CFF-Ball-0x is written in Fortran 95 and requires a fortran compatible compiler. The software comes with a Makefile that compiles the code through the `make` command (`make clean` to clean the compiled files).
The code is parallelised using the `OpenMP` library.
At the moment, the code can be compiled and run on a Unix (or Unix-like) environment.

## How to run the code
The code is simply run by typing in a terminal `./ball0x.exe` (the name of the executable file can be chosen in the Makefile) created post-compilation.
The parameters the code has to be fed with can be found in the self-explanatory file `p.in` (and `param.f90`).
The output file is written in `log.out`, where the current *time unit*, *relative viscosity*, *pressure* and several wall-time parameters are displayed.
The convergence of the simulation, instead, can be checked by plotting the quantities found into `checkConvergence.dat` (first column the *time unit*, other columns some quantities of interest; the second column, for instance, represents the non-normalised relative viscosity. The quantities written in that file can be found inside `timeIntegration.f90`).

## Code acknowledgement
If CFF-Ball-0x has been used, to acknowledge our work please cite the following paper
Monti, A., Rathee, V., Shen, A.Q. and Rosti, M.E., 2021. *A fast and efficient tool to study the rheology of dense suspensions*. Physics of Fluids, 33(10), p.103314. https://doi.org/10.1063/5.0065655

## About us
Alessandro Monti, Ph.D. Postdoctoral Research Fellow at Okinawa Institute of Science and Technology Graduate University (OIST). Creator and main developer
Marco Edoardo Rosti, Ph.D. Assistant Professor at Okinawa Institute of Science and Technology Graduate University (OIST). Coordinator
