export UCX_TLS=tcp
export OMP_NUM_THREADS=1

mkdir -p nest_out/samp

for fgtype in ALG_NE2001 ALG_YMW16
do
    echo "Running nest_samp.py with galaxy type: ${fgtype}"
    mpiexec.hydra -n 4 ./nest_samp.py --fits -g ${fgtype} -o ${fgtype}
done