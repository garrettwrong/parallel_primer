## Building
Note that MPI is a little more involved to build and run.

First you'll actually need an MPI implementation.
A freely available one is OpenMPI, which can usually be installed by

`yum install openmpi openmpi-devel` or equivalent.


MPI implementations are also provided by commercial compiler and computer vendors,
such as Intel and Cray.  At HPC sites, these environments are usually managed by "modules".

In my case, to load openmpi, I run:

`module load mpi`

## Running

Regardless of how you will build, you also need to launch the application using MPI.
There are special commands for this, including many options.
Unfortunately the options names are usually specific to different implementations,
I will be using the OpenMPI flavor.

The most basic is something like the following, where 4 is the number of processes (often called "Ranks").

mpiexec -n 4 ./correlation_dimension_mpi.x

## Remarks

An important feature of MPI is that the processes can be distributed beyond what is possible on a single machine.
Generally speaking, that is what MPI is designed to do, run across many machines in a cluster.
This is in contrast to OpenMP which is targetting a single multicore machine.

That can be critically important for problems that cannot fit into memory on a single machine,
have effecient decompositions, and also for the obvious scaling.
It is generally more effor to write MPI programs,
and can change the "feel" of the code.
This code should probably be refactored to better admit MPI, but then it would be harder to compare the feel implementations.
When considering MPI it is good to consider how both "strong" and "weak" scaling expectations of your problem.

MPI is often combined with OpenMP to form hybrid applications.
