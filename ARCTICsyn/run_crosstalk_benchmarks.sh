#!/bin/bash

NUM_JOBS=2

CIRCDIR=./circuit-structures/complete
LIBDIR=../ARCTICsim/thermo_libs/eval/

LOGFILE=run_$(date +%s%N).txt

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

for f in ${LIBDIR}/*
do
	wait_jobs

	if [[ $f != *.json ]]; then
		continue
	fi

	echo "Processing library $f ..."
	./gradlew -q simulationTestbench --args="-i $CIRCDIR -l $f" </dev/null >> $CIRCDIR/$LOGFILE &
done

for job in `jobs -pr`; do
    wait $job
done

