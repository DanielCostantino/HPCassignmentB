# HPCassignmentB

COMANDI PER COMPILARE:
cd /u/dssc/costa/HPCassignmentB
module load openmpi-4.1.1+gnu-9.3.0
export OMP_PLACES=threads                                                                                                                                                  mpicc -fopenmp density.c -o density.o -lm

RUN:
mpirun --mca pml ucx --mca btl tcp,self --map-by node -np 4 density.o 2 2
