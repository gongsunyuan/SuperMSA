#!/bin/bash
#SBATCH -J msa_opt
#SBATCH -p cnall
#SBATCH -N 1
#SBATCH -o stdout.%j
#SBATCH -e stderr.%j
#SBATCH --no-requeue
#SBATCH --ntasks-per-node=56

# 用法: sbatch <repo>/job_baseline.slurm（stdout/stderr 会被写到提交目录）
# 数据默认取<repo>旁的 ../preliminary，结果写到 <repo>/results；可用环境变量 DATA / OUT / CPU 覆盖。
# 基线约 5–6 小时，分区有时长上限时请加 #SBATCH -t。

SCRIPT=$0
[ -n "${SLURM_JOB_ID:-}" ] && SCRIPT=$(scontrol show job "$SLURM_JOB_ID" | sed -n 's/^ *Command=//p')
ROOT=$(cd "${SLURM_SUBMIT_DIR:-.}" && cd "$(dirname "$SCRIPT")" && pwd)
[ -f "$ROOT/run_queries.sh" ] || { echo "ERROR: 找不到 $ROOT/run_queries.sh" >&2; exit 1; }

# module load compilers/gcc/v12.2.0    # 若编译时用了 module 中的 gcc，在此加载同一环境

# 并发度：默认 7 进程 × 8 线程 = 56 核，可按需覆盖
export JOBS=${JOBS:-7}
export CPU=${CPU:-8}

bash "$ROOT/run_queries.sh"
