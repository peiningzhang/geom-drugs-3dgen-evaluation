#!/bin/bash
export PATH="$PWD/xtb-dist/bin:$PATH"
export PYTHONPATH=/home/phz24002/Drug_Discovery/geom-drugs-3dgen-evaluation:$PYTHONPATH
# 获取时间戳
timestamp=$(date +%Y%m%d_%H%M%S)

# 创建日志文件夹（如果不存在）
mkdir -p logs

# 输入文件路径（默认是 ../semla-flow/elden.sdf，如果传了参数就用参数）
INPUT_SDF=${1:-../semla-flow/elden.sdf}

# 定义日志文件名
LOG_FILE="logs/${timestamp}_energy_pipeline.log"

# 写入日志的函数（带分隔符）
log_section() {
    echo -e "\n==================== $1 ====================\n" >> "$LOG_FILE"
}

# xtb 优化
log_section "xtb_optimization.py"
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 \
python -u scripts/energy_benchmark/xtb_optimization_multi_process.py \
    --input_sdf "$INPUT_SDF" \
    --output_sdf optimized_output.sdf \
    --init_sdf saved_initial_structures.sdf \
    >> "$LOG_FILE" 2>&1

# RMSD 和能量对比
log_section "rmsd_energy.py"
python -u scripts/energy_benchmark/rmsd_energy.py \
    --init_sdf saved_initial_structures.sdf \
    --opt_sdf optimized_output.sdf \
    --n_subsets 5 \
    >> "$LOG_FILE" 2>&1

# 几何对比
log_section "compute_pair_geometry.py"
python -u scripts/energy_benchmark/compute_pair_geometry.py \
    --init_sdf saved_initial_structures.sdf \
    --opt_sdf optimized_output.sdf \
    --n_subsets 5 \
    >> "$LOG_FILE" 2>&1

# 分子稳定性分析（使用 initial 结构）
log_section "compute_molecule_stability.py"
python -u scripts/compute_molecule_stability.py \
    --sdf saved_initial_structures.sdf \
    --aromatic auto \
    --n_subsets 5 \
    >> "$LOG_FILE" 2>&1

echo "All steps completed. Log saved to $LOG_FILE"

