QUARTZ
===================
QUARTZ is a novel and memory-efficient sketch for real-time quantile estimation of persistent items. Intuitively, it identifies candidate persistent items using a bucket-cell multi-level structure, and then employs an associated competitive quantile sketch to profile their value distributions.

## Getting Started
--------------------------------
Tested on Linu.
You need CMake 3.16 or higher and build tools installed.

To compile QUARTZ algorithm:
```
rm -rf build/
mkdir build && cd build
cmake ..
make -j$(nproc)
```
#### Arguments used for testing:
 - **file name** - Path to dataset file.
 - **theta** - Bucket-cell size parameter.
 - **alpha** - Estimate accuracy parameter.
 - **max buckets** - Max bucket size for all DDSketches.
 - **type** -- 1 for memory-THP test;
 -- 2 for persistence filter test;
 -- 3 for tail letency quantile test;
 -- 4 for batched test.
