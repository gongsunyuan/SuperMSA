#!/bin/bash
# 用法: bash build_optimized.sh
# 在测试节点(test01-04)上执行；针对探索1000计算节点的 Xeon Gold 6258R (Cascade Lake) 优化。
# 显式用 -march=cascadelake 而非 native：编译节点与计算节点可能不同，显式指定避免指令集错配。
#
# 构建变体：IMPL=avx512（默认）启用 AVX-512 SSV 滤波器（需要 Xeon Gold 6258R 等
# Cascade Lake 或更新 CPU）；IMPL=sse 回退到原版 SSE 路径。
# LTO 默认关闭：集群系统 ar/ranlib (binutils 2.27) 无法为 gcc12 的 slim LTO 对象
# 建立静态库符号索引，会导致链接时满屏 undefined reference to esl_*。
# 如需尝试 LTO，设 LTO=1：自动追加 -ffat-lto-objects（对象同时含机器码），
# 系统 ar 即可正常建索引，链接阶段仍执行 LTO。
set -euo pipefail

REPO=$(cd "$(dirname "$0")" && pwd)
cd "$REPO"

module load compilers/gcc/v12.2.0 2>/dev/null || true

LTO=${LTO:-0}
IMPL=${IMPL:-avx512}
CFLAGS=${CFLAGS:--O3 -march=cascadelake -funroll-loops -fomit-frame-pointer}
if [ "$LTO" = 1 ]; then
    CFLAGS="$CFLAGS -flto -ffat-lto-objects"
fi

# 先清理旧构建产物，避免被上一次不同参数（尤其 slim LTO）的 .o/.a 污染
if [ -f Makefile ]; then make distclean >/dev/null 2>&1 || true; fi

echo "编译参数: $CFLAGS (IMPL=$IMPL)"
autoconf
if [ "$IMPL" = avx512 ]; then
    ./configure CFLAGS="$CFLAGS" --enable-avx512
else
    ./configure CFLAGS="$CFLAGS"
fi

# 只构建比赛所需目标：easel 静态库 + src（jackhmmer 与 dbprep），跳过 easel/miniapps
make -C easel -j"$(nproc)" libeasel.a
make -C src   -j"$(nproc)"

echo "构建完成，运行冒烟测试..."
./smoke_test.sh
