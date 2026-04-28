#!/bin/bash
# ============================================================
# run_all.sh - 一键运行完整 FRB E_iso 函数分析流程
# 用法: bash run_all.sh [选项]
#
# 选项:
#   cat1      只运行 CHIME Catalog 1 分析
#   cat2      只运行 CHIME Catalog 2 分析
#   repeat    只运行重复爆分析
#   sim       只运行模拟验证
#   plot      只画图
#   (无参数)   运行全部流程
# ============================================================

set -e  # 遇到错误立即停止
export OMP_NUM_THREADS=4

# 颜色输出
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
NC='\033[0m'

log()  { echo -e "${GREEN}[$(date +%H:%M:%S)]${NC} $1"; }
warn() { echo -e "${YELLOW}[$(date +%H:%M:%S)] WARNING:${NC} $1"; }
err()  { echo -e "${RED}[$(date +%H:%M:%S)] ERROR:${NC} $1"; }

# ============================================================
# 0. 环境检查
# ============================================================
log "${CYAN}=== 环境检查 ===${NC}"

# 检查 Python 3
if ! command -v python3 &> /dev/null; then
    err "Python 3 未安装"
    exit 1
fi
PYTHON_VER=$(python3 --version 2>&1)
log "Python: $PYTHON_VER"

# 检查依赖
python3 -c "import numpy, scipy, matplotlib" 2>/dev/null || {
    warn "缺少依赖，正在安装..."
    pip3 install --break-system-packages numpy scipy matplotlib 2>/dev/null || \
    pip3 install numpy scipy matplotlib 2>/dev/null
}

# 检查 PyMultiNest
MPI_AVAILABLE=false
if python3 -c "import pymultinest" 2>/dev/null; then
    if command -v mpiexec &> /dev/null; then
        MPI_AVAILABLE=true
        log "PyMultiNest + MPI: 可用 (完整分析模式)"
    else
        warn "PyMultiNest 已安装但 MPI 不可用"
    fi
else
    warn "PyMultiNest 未安装 (将使用单进程模式，速度较慢)"
fi

# 检查数据文件
for f in chime_cat1.json chime_cat2.json; do
    if [ ! -f "$f" ]; then
        err "缺少数据文件: $f"
        exit 1
    fi
done
log "数据文件: OK"

# 创建输出目录
mkdir -p nest_out/samp nest_out/simu plots/samp plots/simu simu

echo ""

# ============================================================
# 定义运行函数
# ============================================================

run_nest() {
    # 运行 MultiNest 采样
    # $1: 输入文件  $2: 输出前缀  $3: galaxy type  $4+: 额外参数
    local fc="$1" fout="$2" fgt="$3"
    shift 3
    local extra="$@"

    log "运行: nest_samp.py -fc $fc -o $fout -g $fgt $extra"

    if [ "$MPI_AVAILABLE" = true ]; then
        mpiexec -n 4 python3 nest_samp.py -fc "$fc" -o "$fout" -g "$fgt" -upper 1 $extra
    else
        python3 nest_samp.py -fc "$fc" -o "$fout" -g "$fgt" -upper 1 $extra
    fi
}

run_nest_simu() {
    # 运行模拟数据的 MultiNest
    local f1="$1" f2="$2" fout="$3" fgt="$4"

    log "运行: nest_simu.py -f1 $f1 -f2 $f2 -o $fout -g $fgt"

    if [ "$MPI_AVAILABLE" = true ]; then
        mpiexec -n 4 python3 nest_simu.py -f1 "$f1" -f2 "$f2" -o "$fout" -g "$fgt" -upper 1
    else
        python3 nest_simu.py -f1 "$f1" -f2 "$f2" -o "$fout" -g "$fgt" -upper 1
    fi
}

run_plot() {
    # 画后验分布图
    local input="$1" output="$2" title="$3"
    if [ -f "${input}_post_equal_weights.dat" ]; then
        log "画图: $title"
        python3 pltpost.py -f "$input" -o "$output" -title "$title" -up 1 -bo 1
    else
        warn "跳过 $title (结果文件不存在)"
    fi
}

# ============================================================
# 1. CHIME Catalog 1 分析
# ============================================================
run_cat1() {
    log "${CYAN}=== CHIME Catalog 1 分析 (594 FRBs) ===${NC}"

    for fgt in ALG_NE2001 ALG_YMW16; do
        log "Galaxy type: $fgt"
        run_nest "chime_cat1.json" "cat1_${fgt}_upper" "$fgt"
    done

    # 重复爆子样本
    log "分析 Catalog 1 重复爆 (94 bursts)..."
    run_nest "chime_cat1.json" "cat1_repeaters_ALG_YMW16_upper" "ALG_YMW16" "-repeaters"

    # 一次性爆发子样本
    log "分析 Catalog 1 一次性爆发 (500 bursts)..."
    run_nest "chime_cat1.json" "cat1_oneoff_ALG_YMW16_upper" "ALG_YMW16" "-oneoff"

    log "Catalog 1 分析完成"
    echo ""
}

# ============================================================
# 2. CHIME Catalog 2 分析
# ============================================================
run_cat2() {
    log "${CYAN}=== CHIME Catalog 2 分析 (4527 FRBs) ===${NC}"

    for fgt in ALG_NE2001 ALG_YMW16; do
        log "Galaxy type: $fgt"

        # 全部事件
        log "分析全部事件..."
        run_nest "chime_cat2.json" "cat2_${fgt}_upper" "$fgt"

        # 只有一一次性爆发
        log "分析一次性爆发..."
        run_nest "chime_cat2.json" "cat2_oneoff_${fgt}_upper" "$fgt" "-oneoff"
    done

    # 重复爆
    log "分析 Catalog 2 重复爆 (1154 bursts)..."
    run_nest "chime_cat2.json" "cat2_repeaters_ALG_YMW16_upper" "ALG_YMW16" "-repeaters"

    log "Catalog 2 分析完成"
    echo ""
}

# ============================================================
# 3. 模拟验证
# ============================================================
run_sim() {
    log "${CYAN}=== 模拟数据验证 ===${NC}"

    for phis in 1e3 1e4; do
        log "模拟 phis=$phis ..."

        # 生成模拟数据
        python3 simufrb.py -ns 100 -phis ${phis} -alpha -1.5 \
            -logeisomax 41.0 -logeiso0 38.0 \
            -ga 1.4 -bw 400 -ts 25 -fov 200 -sn0 10 \
            -out simu/simdat_${phis}_200.txt

        python3 simufrb.py -ns 100 -phis ${phis} -alpha -1.5 \
            -logeisomax 41.0 -logeiso0 38.0 \
            -ga 0.7 -bw 400 -ts 50 -fov 20 -sn0 10 \
            -out simu/simdat_${phis}_20.txt

        # 恢复参数
        log "恢复 phis=$phis ..."
        run_nest_simu "simu/simdat_${phis}_200.txt" "simu/simdat_${phis}_20.txt" \
            "simdat_${phis}" "ALG_YMW16"

        run_nest_simu "simu/simdat_${phis}_200.txt" "simu/simdat_${phis}_20.txt" \
            "simdat_${phis}_upper" "ALG_YMW16"
    done

    log "模拟验证完成"
    echo ""
}

# ============================================================
# 4. 画图
# ============================================================
run_plots() {
    log "${CYAN}=== 绘制后验分布图 ===${NC}"

    # Catalog 1 结果
    for fgt in ALG_NE2001_upper ALG_YMW16_upper; do
        run_plot "nest_out/samp/cat1_${fgt}" "plots/samp/cat1_${fgt}.eps" \
            "CHIME Cat1 E_iso (${fgt})"
    done
    run_plot "nest_out/samp/cat1_repeaters_ALG_YMW16_upper" \
        "plots/samp/cat1_repeaters.eps" "CHIME Cat1 Repeaters E_iso"
    run_plot "nest_out/samp/cat1_oneoff_ALG_YMW16_upper" \
        "plots/samp/cat1_oneoff.eps" "CHIME Cat1 One-off E_iso"

    # Catalog 2 结果
    for fgt in ALG_NE2001_upper ALG_YMW16_upper; do
        run_plot "nest_out/samp/cat2_${fgt}" "plots/samp/cat2_${fgt}.eps" \
            "CHIME Cat2 E_iso (${fgt})"
        run_plot "nest_out/samp/cat2_oneoff_${fgt}" "plots/samp/cat2_oneoff_${fgt}.eps" \
            "CHIME Cat2 One-off E_iso (${fgt})"
    done
    run_plot "nest_out/samp/cat2_repeaters_ALG_YMW16_upper" \
        "plots/samp/cat2_repeaters.eps" "CHIME Cat2 Repeaters E_iso"

    # 模拟恢复结果
    for phis in 1e3 1e4; do
        run_plot "nest_out/simu/simdat_${phis}" "plots/simu/simdat_${phis}.eps" \
            "Mock data (phis=${phis})"
        run_plot "nest_out/simu/simdat_${phis}_upper" "plots/simu/simdat_${phis}_upper.eps" \
            "Mock data upper (phis=${phis})"
    done

    log "绘图完成"
    echo ""
}

# ============================================================
# 5. 汇总结果
# ============================================================
summary() {
    log "${CYAN}=== 结果汇总 ===${NC}"

    echo ""
    echo "============================================================"
    echo "  FRB E_iso 函数分析结果文件清单"
    echo "============================================================"
    echo ""

    echo "📊 后验采样结果 (nest_out/samp/):"
    ls -lh nest_out/samp/*_stats.dat 2>/dev/null | awk '{print "  " $NF " (" $5 ")"}'
    echo ""

    echo "📈 后验分布图 (plots/samp/):"
    ls -lh plots/samp/*.eps 2>/dev/null | awk '{print "  " $NF " (" $5 ")"}'
    echo ""

    echo "🔬 模拟恢复结果 (nest_out/simu/):"
    ls -lh nest_out/simu/*_stats.dat 2>/dev/null | awk '{print "  " $NF " (" $5 ")"}'
    echo ""

    echo "📊 模拟恢复图 (plots/simu/):"
    ls -lh plots/simu/*.eps 2>/dev/null | awk '{print "  " $NF " (" $5 ")"}'
    echo ""

    # 打印一个 stats 文件的示例
    stats_file=$(ls nest_out/samp/*_stats.dat 2>/dev/null | head -1)
    if [ -n "$stats_file" ]; then
        echo "📋 示例结果 ($stats_file):"
        echo "------------------------------------------------------------"
        cat "$stats_file" 2>/dev/null | head -20
        echo "------------------------------------------------------------"
    fi

    echo ""
    log "${GREEN}全部完成！${NC}"
    echo ""
    echo "关键输出文件:"
    echo "  - 参数统计: nest_out/samp/*_stats.dat"
    echo "  - 后验样本: nest_out/samp/*_post_equal_weights.dat"
    echo "  - 后验图:   plots/samp/*.eps"
    echo ""
}

# ============================================================
# 主流程
# ============================================================
MODE="${1:-all}"

case "$MODE" in
    cat1)
        run_cat1
        run_plot_cat1() {
            for fgt in ALG_NE2001_upper ALG_YMW16_upper; do
                run_plot "nest_out/samp/cat1_${fgt}" "plots/samp/cat1_${fgt}.eps" "CHIME Cat1 E_iso (${fgt})"
            done
            run_plot "nest_out/samp/cat1_repeaters_ALG_YMW16_upper" "plots/samp/cat1_repeaters.eps" "CHIME Cat1 Repeaters"
            run_plot "nest_out/samp/cat1_oneoff_ALG_YMW16_upper" "plots/samp/cat1_oneoff.eps" "CHIME Cat1 One-off"
        }
        run_plot_cat1
        summary
        ;;
    cat2)
        run_cat2
        run_plot_cat2() {
            for fgt in ALG_NE2001_upper ALG_YMW16_upper; do
                run_plot "nest_out/samp/cat2_${fgt}" "plots/samp/cat2_${fgt}.eps" "CHIME Cat2 E_iso (${fgt})"
                run_plot "nest_out/samp/cat2_oneoff_${fgt}" "plots/samp/cat2_oneoff_${fgt}.eps" "CHIME Cat2 One-off (${fgt})"
            done
            run_plot "nest_out/samp/cat2_repeaters_ALG_YMW16_upper" "plots/samp/cat2_repeaters.eps" "CHIME Cat2 Repeaters"
        }
        run_plot_cat2
        summary
        ;;
    repeat)
        log "${CYAN}=== 重复爆专项分析 ===${NC}"
        run_nest "chime_cat1.json" "cat1_repeaters_ALG_YMW16_upper" "ALG_YMW16" "-repeaters"
        run_nest "chime_cat2.json" "cat2_repeaters_ALG_YMW16_upper" "ALG_YMW16" "-repeaters"
        run_plot "nest_out/samp/cat1_repeaters_ALG_YMW16_upper" "plots/samp/cat1_repeaters.eps" "CHIME Cat1 Repeaters"
        run_plot "nest_out/samp/cat2_repeaters_ALG_YMW16_upper" "plots/samp/cat2_repeaters.eps" "CHIME Cat2 Repeaters"
        summary
        ;;
    sim)
        run_sim
        run_plots_sim() {
            for phis in 1e3 1e4; do
                run_plot "nest_out/simu/simdat_${phis}" "plots/simu/simdat_${phis}.eps" "Mock (phis=${phis})"
                run_plot "nest_out/simu/simdat_${phis}_upper" "plots/simu/simdat_${phis}_upper.eps" "Mock upper (phis=${phis})"
            done
        }
        run_plots_sim
        summary
        ;;
    plot)
        run_plots
        summary
        ;;
    all)
        log "${CYAN}============================================${NC}"
        log "${CYAN}  FRB E_iso 函数完整分析流程${NC}"
        log "${CYAN}  数据: CHIME Catalog 1 + 2 + Repeaters${NC}"
        log "${CYAN}============================================${NC}"
        echo ""
        START_TIME=$(date +%s)

        run_cat1
        run_cat2
        run_sim
        run_plots

        END_TIME=$(date +%s)
        ELAPSED=$((END_TIME - START_TIME))
        log "总耗时: $((ELAPSED/60))分$((ELAPSED%60))秒"

        summary
        ;;
    *)
        echo "用法: bash run_all.sh [cat1|cat2|repeat|sim|plot|all]"
        echo ""
        echo "选项:"
        echo "  cat1     只运行 CHIME Catalog 1 分析"
        echo "  cat2     只运行 CHIME Catalog 2 分析"
        echo "  repeat   只运行重复爆分析"
        echo "  sim      只运行模拟验证"
        echo "  plot     只画图"
        echo "  all      运行全部流程 (默认)"
        exit 1
        ;;
esac
