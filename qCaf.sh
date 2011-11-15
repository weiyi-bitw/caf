#!/bin/sh
#$ -j y
#$ -cwd
#$ -S /bin/bash
#$ -t 1-700
export CONF=$1
cd /home/weiyi/javaworks/caf
/usr/java/jdk1.6.0_16/bin/java -Xmx2000M -DSGE_TASK_LAST=$SGE_TASK_LAST -DSGE_TASK_ID=$SGE_TASK_ID -DJOB_ID=$JOB_ID -cp ./dist/caf-0.2.jar caf.CorrAttractorFinder $CONF
