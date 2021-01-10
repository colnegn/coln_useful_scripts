#!/bin/bash


## HARD-CODED path to test output file:
test_out='output'
#test_out='a__o2h3/output'


if [ $# -ne 3 ]; then
  echo "usage: 1.input (FULL PATH) directory list;  2.output directory list for jobs to resubmit; 3.test-output path from dlist directories"
  exit
fi

in_dlist=$1
out_dlist=$2
test_out=$3



if [ -f $out_dlist ]; then
  echo "ERROR: $out_dlist already exists... exiting..."
  exit
fi


for d in $(cat $in_dlist); do

  ## check if $test_out exists:
  if ! [ -f $d/$test_out ]; then
    ## no output file, so add to output directory list:
    echo "$d" >> $out_dlist
    continue
  fi
  
  ## check if job (represented by $test_out calculation) completed successfully:
  jobdone=`grep -c 'Thank you very much for using Q-Chem' $d/$test_out`
  
  ## calculation completed:
  if [ $jobdone -ne 0 ]; then
    ## print job directory to successful dlist:
    echo "$d" >> ${out_dlist}_SUCCESS
    ## iterate for-loop:
    continue
  fi
    
  ## if calculation did not complete, see if it passed license check:
  licensecheck=`grep -c 'SCF failed to converge' $d/$test_out`

  ## failed license check:
  if [ $licensecheck -eq 0 ]; then
    ## so print job directory path to output directory list:
    echo "$d" >> $out_dlist

  ## else job just failed normally (?):
  else
    ## print job directory to failed dlist:
    echo "$d" >> ${out_dlist}_FAILED
  fi

done













