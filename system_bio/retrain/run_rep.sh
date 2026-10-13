#!/bin/bash
cd "D:/Cander/School/3-2/system_bio"
for job in "ss 36" "aug 36" "ctl1 12" "ctl2 12" "ctl3 12"; do
  set -- $job
  if [ -f retrain/out/balanced_$1_$2.pkl ]; then echo "skip $1 $2"; continue; fi
  tools/pdl/Scripts/python.exe retrain/pipeline.py balanced $1 $2 2 2>&1 | grep -v -i "executable"
done
echo REP_DONE
