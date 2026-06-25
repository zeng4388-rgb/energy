export UCX_TLS=tcp
export OMP_NUM_THREADS=1

fov1=0.55
fov2=30
galaxy_type=ALG_YMW16

mkdir -p nest_out/simu
mkdir -p simu

for phis in 1e3 1e4
do
    fout1=simdat_${phis}
    fout2=simdat_${phis}_upper
    fin1=./simu/simdat_${phis}_${fov1}.txt
    fin2=./simu/simdat_${phis}_${fov2}.txt

    # 检查输入文件是否存在
    if [ ! -f "$fin1" ] || [ ! -f "$fin2" ]; then
        echo "ERROR: 缺少模拟输入文件，请先运行 simudat.sh 生成数据"
        echo "  缺失: $fin1 或 $fin2"
        exit 1
    fi

    echo "Running nest_simu.py for phi*=${phis}"
    mpiexec.hydra -n 4 ./nest_simu.py -f1 $fin1 -f2 $fin2 -o $fout1 -g ${galaxy_type}
    mpiexec.hydra -n 4 ./nest_simu.py -f1 $fin1 -f2 $fin2 -o $fout2 -g ${galaxy_type}
done
exit
