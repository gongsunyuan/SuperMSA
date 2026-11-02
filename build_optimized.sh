#!/bin/bash
# 用法: bash build_optimized.sh
# 在测试节点(test01-04)上执行；针对探索1000计算节点的 Xeon Gold 6258R (Cascade Lake) 优化。
# 显式用 -march=cascadelake 而非 native：编译节点与计算节点可能不同，显式指定避免指令集错配。
# 若 LTO 链接失败，设 LTO=0 重试；若需对照基线，设 CFLAGS 环境变量覆盖。
set -euo pipefail

REPO=$(cd "$(dirname "$0")" && pwd)
cd "$REPO"

module load compilers/gcc/v12.2.0 2>/dev/null || true

LTO=${LTO:-1}
CFLAGS=${CFLAGS:--O3 -march=cascadelake -funroll-loops -fomit-frame-pointer}
if [ "$LTO" = 1 ]; then
    CFLAGS="$CFLAGS -flto"
fi

echo "编译参数: $CFLAGS"
autoconf
./configure CFLAGS="$CFLAGS"
make -j"$(nproc)"

echo "构建完成，运行冒烟测试..."
./smoke_test.sh
