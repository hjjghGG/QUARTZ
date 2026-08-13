QUARTZ
===================
QUARTZ is a novel and memory-efficient sketch for real-time quantile estimation of persistent items. Intuitively, it identifies candidate persistent items using a bucket-cell multi-level structure, and then employs an associated competitive quantile sketch to profile their value distributions.

## Getting Started
--------------------------------
Tested on Linux.
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

#### module enable/disable switch and parameter tuning
 - **RS_IN_CELL** - enable/disable Reservoir Sample in Cell optimization.
 - **ORG_COLLAPSE** - enable/disable original collapse strategy instead of the adaptive collapse strategy.
 - **KCU_CELL_SIZE** - number of cells in a bucket (excluding the champion).
 - **KCU_EPOCH_SIZE** - time window size.
 - **KCU_epsilon** -- epsilon.
 - **KCU_epsilon2** -- miu.
 - Replacing the KC_pu file with KC_pu_improved and modifying the relevant CMakeLists.txt enables Shared Persistence Filter.
 
