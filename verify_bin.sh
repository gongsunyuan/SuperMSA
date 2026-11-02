#!/bin/bash
# 用法: bash verify_bin.sh   （在测试节点 test01-04 上运行）
# 验证 PMSA 二进制数据库与 FASTA 产生完全一致的 jackhmmer 结果。
# 默认用 data/demodb.fasta + data/query.fasta；可用 DB/QUERY/CPU 环境变量覆盖，
# 建议也抽一个正式 query + 正式库跑一次（DB=../preliminary/database.fasta QUERY=...）。
set -euo pipefail

REPO=$(cd "$(dirname "$0")" && pwd)
DB=${DB:-$REPO/data/demodb.fasta}
QUERY=${QUERY:-$REPO/data/query.fasta}
JACKHMMER=$REPO/src/jackhmmer
DBPREP=$REPO/src/dbprep
CPU=${CPU:-8}
FILTERS="--noali --F1 0.0005 --F2 5e-05 --F3 5e-07"

[ -x "$JACKHMMER" ] || { echo "ERROR: 先编译 jackhmmer" >&2; exit 1; }
[ -x "$DBPREP" ]    || { echo "ERROR: 先编译 dbprep" >&2; exit 1; }

WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

"$DBPREP" "$DB" "$WORK/db.bin"

for fmt in fasta bin; do
    if [ "$fmt" = bin ]; then dbfile="$WORK/db.bin"; else dbfile="$DB"; fi
    "$JACKHMMER" --cpu "$CPU" $FILTERS -N 3 -E 0.001 \
        --tblout "$WORK/$fmt.tbl" -A "$WORK/$fmt.sto" -o "$WORK/$fmt.out" \
        "$QUERY" "$dbfile"
done

if diff <(grep -v '^#' "$WORK/fasta.tbl") <(grep -v '^#' "$WORK/bin.tbl"); then
    echo "PASS: tbl 数据行完全一致"
else
    echo "FAIL: tbl 不一致" >&2; exit 1
fi

if diff "$WORK/fasta.sto" "$WORK/bin.sto" > /dev/null; then
    echo "PASS: msa (sto) 完全一致"
else
    echo "FAIL: sto 不一致" >&2; exit 1
fi

echo "全部通过：二进制库结果与 FASTA 逐字节一致"
