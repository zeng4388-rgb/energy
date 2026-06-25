#!/bin/sh

# 绘制真实样本后验分布图
# 输出目录: plots/samp/

mkdir -p plots/samp

for fgtype in ALG_NE2001 ALG_YMW16
do
    echo "Plotting posterior for ${fgtype}..."
    python3 pltpost.py -f ./nest_out/samp/${fgtype} -o ./plots/samp/${fgtype}.eps -title "Real Sample - ${fgtype}" -up -bo
done
