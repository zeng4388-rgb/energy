#!/bin/sh

# 绘制模拟数据后验分布图
# 输出目录: plots/simu/

mkdir -p plots/simu

for phis in 1e3 1e4
do
    o1=./plots/simu/simdat_${phis}.eps
    o2=./plots/simu/simdat_${phis}_upper.eps
    in1=./nest_out/simu/simdat_${phis}
    in2=./nest_out/simu/simdat_${phis}_upper
    python3 pltpost.py -f $in1 -o $o1 -title "Mock data phi*=${phis}" -up -bo
    python3 pltpost.py -f $in2 -o $o2 -title "Mock data phi*=${phis} upper" -up -bo
done
exit
