#!/bin/bash
cd "D:/Cander/School/3-2/system_bio"
for mt in balanced topology mixed; do
  for kind in ss rand; do
    if [ -f retrain/out/${mt}_${kind}_24.pkl ]; then echo "skip $mt $kind (done)"; continue; fi
    tools/pdl/Scripts/python.exe retrain/pipeline.py $mt $kind 24 5 2>&1 | grep -v -i "executable"
  done
done
echo BIG_DONE
