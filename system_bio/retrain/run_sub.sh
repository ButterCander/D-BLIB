#!/bin/bash
cd "D:/Cander/School/3-2/system_bio"
for kind in sub12_1 sub12_2; do
  if [ -f retrain/out/balanced_${kind}_12.pkl ]; then echo "skip $kind"; continue; fi
  tools/pdl/Scripts/python.exe retrain/pipeline.py balanced $kind 12 6 2>&1 | grep -v -i "executable"
done
echo SUB_DONE
