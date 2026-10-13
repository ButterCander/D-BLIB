#!/bin/bash
cd "D:/Cander/School/3-2/system_bio"
for r in 12 24; do
  if [ -f retrain/out/balanced_aug_${r}.pkl ]; then echo "skip aug $r"; continue; fi
  tools/pdl/Scripts/python.exe retrain/pipeline.py balanced aug $r 6 2>&1 | grep -v -i "executable"
done
echo AUG_DONE
