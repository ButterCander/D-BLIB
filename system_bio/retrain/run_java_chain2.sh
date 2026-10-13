#!/bin/bash
cd "D:/Cander/School/3-2/system_bio"
export JXMX=1000m
bash retrain/run_java.sh tl_ss ss 3 10 30
bash retrain/run_java.sh tl_rand rand 3 10 30
echo JAVA_CHAIN2_DONE
