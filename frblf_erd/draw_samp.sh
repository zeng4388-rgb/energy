#!/usr/bin/env bash

# 切换到工作目录
cd "$(dirname "$0")"

# 确保输出目录存在
mkdir -p ./plots

for fgtype in ALG_NE2001_upper ALG_YMW16_upper
do
    echo "New plot of ${fgtype} has been finished."
    python3 ./pltpost.py -f ./nest_out/samp/${fgtype} -o ./plots/${fgtype}.eps -title "Real Sample" -up 1 -bo 1
done
