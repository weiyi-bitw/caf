#!/bin/sh
#$ -j y
#$ -cwd
#$ -S /bin/bash
#$ -t 1-600
export CONF=$1
export COMMAND=$2
export JOBID=$3 # JOBID for debugging mode / attractor folder for MRC mode
export BREAKPOINT=$4 # break point for debugging mode / minsize for MRC mode
export OVLPTH=$5 # overlap threshold for MRC mode
cd /home/weiyi/javaworks/caf
/usr/java/jdk1.6.0_16/bin/java -Xmx2000M -DSGE_TASK_LAST=$SGE_TASK_LAST -DSGE_TASK_ID=$SGE_TASK_ID -DJOB_ID=$JOB_ID -cp ./dist/caf-0.2.jar caf.CorrAttractorFinder $CONF $COMMAND $JOBID $BREAKPOINT $OVLPTH
