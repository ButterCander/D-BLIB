#!/bin/bash
# step 1 in resumable chunks of 12 GA runs (36 models): seeds 100-111, 112-123, 124-135 (disjoint from all earlier batches); thread backend, 3 threads
cd "D:/Cander/School/3-2/system_bio"
for base in 124 112 100; do
  for kind in rand ss; do
    f=retrain/out/balanced_${kind}_12_b${base}.pkl
    if [ -f $f ]; then echo "skip $kind $base"; continue; fi
    THREADS=1 SEED_BASE=$base tools/pdl/Scripts/python.exe retrain/pipeline.py balanced $kind 12 3 2>&1 | grep -v -i "executable"
  done
done
echo BIG2B_DONE
