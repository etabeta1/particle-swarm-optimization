#PBS -q cpu
#PBS -N swarm_search_benchmark
#PBS -l select=1:ncpus=28
#PBS -l walltime=08:00:00
#PBS -j oe

cd $PBS_O_WORKDIR
make clean && make dist
find ./bin -type f -exec ./{} \;
