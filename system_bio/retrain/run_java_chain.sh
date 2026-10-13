#!/bin/bash
cd "D:/Cander/School/3-2/system_bio"
export JXMX=1200m
bash retrain/run_java.sh topo_ss ss 0 10 30
bash retrain/run_java.sh topo_rand rand 0 10 30
echo JAVA_CHAIN_DONE
