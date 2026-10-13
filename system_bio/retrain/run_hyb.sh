#!/bin/bash
cd "D:/Cander/School/3-2/system_bio"
for mt in balanced topology mixed; do
  if [ -f retrain/out/${mt}_hyb_12.pkl ]; then echo "skip $mt hyb (done)"; continue; fi
  tools/pdl/Scripts/python.exe retrain/pipeline.py $mt hyb 12 6 2>&1 | grep -v -i "executable"
done
echo HYB_DONE
