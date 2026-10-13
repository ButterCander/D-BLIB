#!/bin/bash
# usage: run_java.sh <name> <ss|rand> <balance> <topology> <simulations>   (Java Gitsbe+Drabme pipeline, CASCADE 2.0)
cd "D:/Cander/School/3-2/system_bio"
name=$1; kind=$2; bal=$3; topo=$4; sims=$5
dir=retrain/java_runs/$name
rm -rf $dir && mkdir -p $dir && cp tools/java/druglogics-synergy/ags_cascade_2.0/{config,drugpanel,modeloutputs,network.sif,perturbations} $dir/
if [ "$kind" = "ss" ]; then cp tools/java/druglogics-synergy/ags_cascade_2.0/steadystate $dir/training; else cp tools/java/druglogics-synergy/ags_cascade_2.0/random_train $dir/training; fi
sed -i "s/^simulations:.*/simulations:\t$sims/; s/^balance_mutations:.*/balance_mutations:\t$bal/; s/^topology_mutations:.*/topology_mutations:\t$topo/; s/^synergy_method:.*/synergy_method:\tbliss/; s/^parallel_sim_num:.*/parallel_sim_num:\t4/; s/^verbosity:.*/verbosity:\t0/" $dir/config
export JAVA_HOME="C:/Program Files/Microsoft/jdk-17.0.20.101-hotspot"
start=$(date +%s)
"$JAVA_HOME/bin/java" -Xmx${JXMX:-1400m} -cp tools/java/druglogics-synergy/target/synergy-1.2.1-jar-with-dependencies.jar eu.druglogics.synergy.Launcher --inputDir=$dir --project=$name > $dir/run.log 2>&1
echo "$name ($kind, balance=$bal topology=$topo sims=$sims): exit $?, $(( $(date +%s) - start )) s"
