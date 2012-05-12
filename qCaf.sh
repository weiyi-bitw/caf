#!/bin/sh
#$ -j y
#$ -cwd
#$ -S /bin/bash
#$ -t 1-7000
export CONF=$1
export JOBID=$2
export BREAKPOINT=$3
cd /home/weiyi/workspace/javaworks/caf
java -Xmx2000M -DSGE_TASK_LAST=7000 -DSGE_TASK_ID=460 -DJOB_ID=101 -cp ./dist/caf-0.2.jar caf.CorrAttractorFinder $CONF $JOBID $BREAKPOINT
