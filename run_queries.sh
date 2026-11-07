#!/bin/bash
# 用法: run_queries.sh [data_dir] [out_dir]
#   data_dir  含 query/*.fasta 与 database.fasta，默认 ../preliminary（与仓库并列）
#   out_dir   默认 <repo>/results，也可用环境变量 OUT 指定
# 输出:
#   tbl/<sid>.tbl   命中表（评分输入）     msa/<sid>.sto  多序列比对
#   timing.tsv      逐条耗时（仅供参考）   wall_time.tsv  正式计时：预热结束 → 最后一条跑完
# 环境变量: JACKHMMER CPU(8) JOBS(56/CPU) ITERS(3) EVALUE(0.001) FILTERS WARMUP(1) FORCE(0，设 1 重跑已有结果)
#          NUMA(0，设 1 内存交织) AFFINITY(0，设 socket 按双路绑定) NUMA_NODES(2)
#          ORDER(longest，name 恢复文件名顺序) HISTORY(可选历史 timing.tsv，默认 results_sse/timing.tsv)
# 并发: 同时运行 JOBS 个 jackhmmer 进程（每个占 CPU 线程），动态队列，跑完一条立即补下一条。
#       默认 7 进程 × 8 线程 = 56 核，与探索1000 单节点匹配。
# 预热: 先对第一条 query 跑 WARMUP 次，把数据库读进 page cache，结果丢弃、不计时。
#       一次完整调用已含多轮全库扫描，足以灌满 page cache，故默认 1。
# NUMA: 设 1 时各进程以 numactl --interleave=all 启动（双 socket 节点内存交织）。
# 评分: data_dir 下有 baseline_roc10.tsv 时，逐条判定 ROC10 是否不低于基线的 90%。
set -euo pipefail

case $(hostname) in
    cln0[1-4]*|ibcln*)
        echo "ERROR: 禁止在登录节点运行程序。小规模测试请 ssh test01；完整运行请用 sbatch 提交（见 job_baseline.slurm）。" >&2
        exit 1;;
esac

REPO=$(cd "$(dirname "$0")" && pwd)
DATA=${1:-${DATA:-$(dirname "$REPO")/preliminary}}
OUT=${2:-${OUT:-$REPO/results}}
JACKHMMER=${JACKHMMER:-$REPO/src/jackhmmer}
CPU=${CPU:-8}
JOBS=${JOBS:-$((56 / CPU))}
ITERS=${ITERS:-3}
EVALUE=${EVALUE:-0.001}
FILTERS=${FILTERS:---noali --F1 0.0005 --F2 5e-05 --F3 5e-07}
WARMUP=${WARMUP:-1}
FORCE=${FORCE:-0}
NUMA=${NUMA:-0}
AFFINITY=${AFFINITY:-0}
NUMA_NODES=${NUMA_NODES:-2}
ORDER=${ORDER:-longest}
HISTORY=${HISTORY:-$REPO/results_sse/timing.tsv}
DB=$DATA/database.fasta
TIMING=$OUT/timing.tsv

die() { echo "ERROR: $*" >&2; exit 1; }

[ "$JOBS" -ge 1 ] || die "JOBS 必须 >= 1（当前 CPU=$CPU，请显式设 JOBS）"
[ "$ORDER" = longest ] || [ "$ORDER" = name ] || die "ORDER 仅支持 longest 或 name"
if [ "$AFFINITY" != 0 ]; then
    [ "$AFFINITY" = socket ] || die "AFFINITY 仅支持 0 或 socket"
    [ "$NUMA" = 0 ] || die "AFFINITY=socket 与 NUMA=1 不可同时使用"
    [ "$NUMA_NODES" -ge 1 ] || die "NUMA_NODES 必须 >= 1"
    command -v numactl >/dev/null || die "AFFINITY=socket 需要 numactl"
    xargs --help 2>&1 | grep -q -- --process-slot-var || die "当前 xargs 不支持 --process-slot-var"
fi
if [ $((JOBS * CPU)) -gt 56 ]; then
    echo "WARNING: JOBS*CPU=$((JOBS * CPU)) 超过 56 核，将超订 CPU" >&2
fi

[ -f "$DB" ]           || die "$DB not found"
[ -d "$DATA/query" ]   || die "$DATA/query not found"
[ -x "$JACKHMMER" ]    || die "jackhmmer not executable: $JACKHMMER"
queries=("$DATA"/query/*.fasta)
[ -e "${queries[0]}" ] || die "no *.fasta in $DATA/query"

mkdir -p "$OUT/tbl" "$OUT/msa" "$OUT/timing" 2>/dev/null && touch "$OUT/.write_test" 2>/dev/null \
    || die "输出目录不可写: $OUT（home 在计算节点只读，请把仓库放在 ~/WORK 下或指定 out_dir）"
rm -f "$OUT/.write_test"
if [ "$FORCE" = 1 ]; then rm -f "$OUT"/timing/*.time; fi

# 数据库预处理：一次性转换为 PMSABIN1 二进制格式（mmap 读取，跳过 FASTA 解析）。
# 仅在 .bin 缺失或比 fasta 旧时转换；耗时追加到 $OUT/preprocess.tsv（正式报告需公开此时间）。
DBPREP=$REPO/src/dbprep
if [ -x "$DBPREP" ]; then
    BINDB=$DATA/database.bin
    if [ ! -f "$BINDB" ] || [ "$DB" -nt "$BINDB" ]; then
        echo "预处理: dbprep $DB -> $BINDB"
        t_prep=$(date +%s.%N)
        "$DBPREP" "$DB" "$BINDB"
        prep_secs=$(awk -v a="$t_prep" -v b="$(date +%s.%N)" 'BEGIN{printf "%.1f", b-a}')
        printf 'preprocess\t%s\n' "$prep_secs" >> "$OUT/preprocess.tsv"
        echo "预处理耗时: ${prep_secs} s（一次性，摊薄到全部 query）"
    fi
    DB=$BINDB
fi

ARGS=(--cpu "$CPU" $FILTERS -N "$ITERS" -E "$EVALUE" -o /dev/null)
export JACKHMMER DB OUT NUMA AFFINITY NUMA_NODES
export ARGS_STR="${ARGS[*]}"

run_one() {
    local qf=$1
    local sid secs hits t0 t1
    sid=$(basename "$qf" .fasta)
    t0=$(date +%s.%N)
    if [ "$AFFINITY" = socket ]; then
        node=$(( ${PMSA_SLOT:-0} % NUMA_NODES ))
        numactl --cpunodebind="$node" --membind="$node" "$JACKHMMER" $ARGS_STR \
            --tblout "$OUT/tbl/$sid.tbl" -A "$OUT/msa/$sid.sto" "$qf" "$DB"
    elif [ "$NUMA" = 1 ]; then
        numactl --interleave=all "$JACKHMMER" $ARGS_STR \
            --tblout "$OUT/tbl/$sid.tbl" -A "$OUT/msa/$sid.sto" "$qf" "$DB"
    else
        "$JACKHMMER" $ARGS_STR \
            --tblout "$OUT/tbl/$sid.tbl" -A "$OUT/msa/$sid.sto" "$qf" "$DB"
    fi
    t1=$(date +%s.%N)
    secs=$(awk -v a="$t0" -v b="$t1" 'BEGIN{printf "%.1f", b-a}')
    hits=$(grep -vc '^#' "$OUT/tbl/$sid.tbl" || true)
    printf '%s\t%s\t%s\n' "$sid" "$secs" "$hits" > "$OUT/timing/$sid.time"
    echo "$sid: ${secs}s, $hits hits"
}
export -f run_one

total=${#queries[@]}

if [ "$WARMUP" -gt 0 ]; then
    echo "预热 $WARMUP 次（不计时）..."
    for _ in $(seq "$WARMUP"); do
        "$JACKHMMER" "${ARGS[@]}" --tblout /dev/null "${queries[0]}" "$DB"
    done
fi

# 筛选待办：已有非空 tbl 且 FORCE=0 的跳过
TODO=$(mktemp "$OUT/.todo.XXXXXX")
trap 'rm -f "$TODO"' EXIT
skipped=0
for qf in "${queries[@]}"; do
    sid=$(basename "$qf" .fasta)
    if [ "$FORCE" = 0 ] && [ -s "$OUT/tbl/$sid.tbl" ]; then
        skipped=$((skipped + 1))
        continue
    fi
    printf '%s\n' "$qf" >> "$TODO"
done
ran=$((total - skipped))
echo "共 $total 条 query：本次运行 $ran 条，跳过 $skipped 条；并发 $JOBS 进程 × $CPU 线程"

if [ "$ORDER" = longest ]; then
    HISTORY_FILE=
    [ ! -f "$HISTORY" ] || HISTORY_FILE=$HISTORY
    ORDERED=$(mktemp "$OUT/.ordered.XXXXXX")
    trap 'rm -f "$TODO" "$ORDERED"' EXIT
    awk -v hist="$HISTORY_FILE" -v list="$TODO" '
        BEGIN {
            while (hist != "" && (getline line < hist) > 0) {
                n = split(line, field, "\t")
                if (n >= 2 && field[2] + 0 > 0) seconds[field[1]] = field[2] + 0
            }
            while ((getline qf < list) > 0) {
                id = ""; len = 0
                while ((getline line < qf) > 0) {
                    if (line ~ /^>/) {
                        if (id == "") {
                            id = substr(line, 2)
                            sub(/[[:space:]].*$/, "", id)
                        }
                    } else {
                        gsub(/[[:space:]]/, "", line)
                        len += length(line)
                    }
                }
                close(qf)
                if (id != "") printf "%.3f\t%d\t%s\n", seconds[id] + 0, len, qf
            }
        }
    ' /dev/null | sort -t $'\t' -k1,1nr -k2,2nr -k3,3 | cut -f3 > "$ORDERED"
    mv "$ORDERED" "$TODO"
fi

echo "队列顺序: $ORDER$([ "$ORDER" = longest ] && [ -n "${HISTORY_FILE:-}" ] && echo "（历史: $HISTORY_FILE）" || true)"
t_all=$(date +%s.%N)
if ! xargs -a "$TODO" -r -d '\n' -n 1 -P "$JOBS" --process-slot-var=PMSA_SLOT bash -c 'run_one "$@"' _; then
    die "有 query 运行失败，请检查上方日志"
fi
wall=$(awk -v a="$t_all" -v b="$(date +%s.%N)" 'BEGIN{printf "%.1f", b-a}')

# 汇总逐条耗时（含以往运行留下的记录，与"仅供参考"定位一致）
{
    echo -e "sid\tseconds\thits"
    cat "$OUT"/timing/*.time 2>/dev/null || true
} > "$TIMING"

echo -e "seconds\tqueries_run\tskipped\n$wall\t$ran\t$skipped" > "$OUT/wall_time.tsv"
echo "正式计时: ${wall} s（运行 $ran 条，跳过 $skipped 条）"
[ "$skipped" = 0 ] || echo "WARNING: 有 $skipped 条沿用已有结果，计时无效；正式测试请清空输出目录或设 FORCE=1" >&2

BASELINE=$DATA/baseline_roc10.tsv
[ -f "$BASELINE" ] || BASELINE=
python3 "$REPO/score.py" --data "$DATA" --results "$OUT" ${BASELINE:+--baseline "$BASELINE"} || echo "warning: scoring step failed" >&2
