# PhD
## MPI Fortran code to solve Poisson's equation using Gauss-Seidel algorithm (or Red-Black algorithm)

- 'Makefile': used to compile the Fortran file and generate the executable. (MPI module should be loaded first in the parallel case)
- 'serial.f90': first version of the code in serial, without MPI.
- 'parallel.f90': parallel version of the code to solve Poisson's equation for a given (punctual or gaussian) charge distribution and a boundary conditions defined. It writes the source, the potential obtained with Gauss-Seidel (GS) and the potential obtained with Red-Black (RB) in separate files.
- 'plot-src.py': simple Python code to plot the source (right hand side in Poisson's equation) by reading the file 'src.dat' generated with the Fortran executable.
- 'plot-gs.py': simple Python code to plot the potential obtained with GS algorithm by reading the file 'pot_gs.dat' generated with the Fortran executable.
- 'plot-rb.py': simple Python code to plot the potential obtained with RB algorithm by reading the file 'pot_rb.dat' generated with the Fortran executable.

Example of how run the executable file generated with the Fortran compilation with 8 tasks: mpirun -np 8 ./parallel > job.out
(in this case, 'job.out' will contain the number of iterations and time needed by both algorithms to solve Poisson's equation)

More details about the code are available as comments through the 'parallel.f90' file.