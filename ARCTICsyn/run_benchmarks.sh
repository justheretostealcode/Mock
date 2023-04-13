#!/bin/bash

NUM_JOBS=1

OUTDIR=benchmarks
LOGFILE=enumeration_$(date +%s%N).txt

function handler {
    kill -- -$$
}

trap "handler" SIGINT

function wait_jobs {
  RUNNING_JOBS=($(jobs -pr))
	RUNNING_JOBS=${#RUNNING_JOBS[*]}
	while [[ $RUNNING_JOBS -ge ${NUM_JOBS} ]]; do
		sleep 1
		RUNNING_JOBS=($(jobs -pr))
		RUNNING_JOBS=${#RUNNING_JOBS[*]}
	done
}

export GRADLE_OPTS=-Xmx20G

while IFS=\n read -r truthtable
do
	wait_jobs
    echo "mapping truth table $truthtable ..."
	./gradlew -q run --args="-t $truthtable" </dev/null >> $OUTDIR/$LOGFILE &
done < cello_functions.csv

for job in `jobs -pr`; do
    wait $job
done
