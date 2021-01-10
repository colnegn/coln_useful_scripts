#!/bin/bash

#!# time to sleep while waiting for calculations:
sleeptime='30m'

#!# do we have MO and Hessian files already?
MO_initrestart_bool=1
HESS_initrestart_bool=0

## print out hard-coded warning:
echo "HARD-CODED: sleep time = $sleeptime ;  MO restart bool = $MO_initrestart_bool ; HESS restart bool = $HESS_initrestart_bool ... sleeping for 3s to allow you to kill job..."
sleep 5s


if [ $# -ne 5 ]; then
  echo "usage: 1.input file (with molecule, scf_guess, and hessian = read); 2.scratch directory (containing molecule, and other scratch files if restart); 3.directory prefix (for directories to be created); 4.submission script; 5.output file name (in submission script)"
  exit
fi

initinput="$1"
scratchdir="$2"
dirprefix="$3"
subscript="$4"
outfilename="$5"



if ! [ -f ./write_fake_hessian.py ]; then
  echo "ERROR: write_fake_hessian.py must be in working directory... exiting..."
  exit
fi

## make sure we don't overwrite directories:
if [ $(ls | grep -c "$dirprefix") -gt 0 ]; then
  echo "ERROR: remove files/directories with name $dirprefix ... exiting..."
  exit
fi

module load python


here=`pwd`

## calculation counter:
count=1
## setup initial optimization:
d="${dirprefix}_${count}"
mkdir "$d"
cp -r "$initinput" "$scratchdir" "$subscript" "$d"
#cd "$d"


## do we already have MOs?
if [ $MO_initrestart_bool != 1 ]; then
  ## don't have orbitals or Hessian to read:
  sed -i 's/scf_guess/!scf_guess/' "$d/$initinput"
fi

## do we already have HESS (BFGS guess)?
if [ $HESS_initrestart_bool != 1 ]; then
  ## don't have orbitals or Hessian to read:
  sed -i 's/GEOM_OPT_HESSIAN/!GEOM_OPT_HESSIAN/' "$d/$initinput"
fi


#!# make temp dir:
#cd "$here"
cp -r "$d" ./tmp_${d}
cd "$d"

## submit initial job:
qsub -N "${count}_${dirprefix}" "$subscript" > job.id

cd "$here"


jobscomplete=0
## start main loop:
while [ $jobscomplete != 1 ]; do

  ## look for file which is only written after calculation completes:
  if [ -f "$d/STDOUT" ]; then
    
    ## check if optimization converged:
    optconverged_str='Thank you very much for using Q-Chem'
    if [ $(grep -c "$optconverged_str" "$d/$outfilename") -eq 1 ]; then
      ## holy shit it actually worked????
      echo "optimization converged?"
      jobscomplete=1
      exit
    fi


    #!# check if we DO NOT have *NEW* MO guesses:
    if ! [ "$(diff $d/$scratchdir/53.0 ./tmp_${d}/$scratchdir/53.0 | wc -l)" -eq 0 ]; then
    ##EVALS TRUE IF 53.0 FROM PREV, BUT NO ITERS:if ! [ -f "$d/$scratchdir/53.0" ]; then
      #!# keep recording of errors:
      echo "failure in $d"
      cat $d/STD*
      tail -200 "$d/$outfilename"
      echo

      ## restart previous calculation:
      rm -r "$d"
      cp -r ./tmp_${d} "$d"
      cd "$d"
      qsub -N "${count}_${dirprefix}" "$subscript" > job.id
      cd "$here"
      ## back to while loop:
      continue
    fi

    #!# if we made it here, do next calculation:
    dprev="$d"
    count=$((count + 1))
    d="${dirprefix}_${count}"
    mkdir "$d"
    cp "$initinput" "$subscript" "$d"
    mkdir "$d/$scratchdir"
    cp "$dprev/$scratchdir/molecule" "$d/$scratchdir"
    cp "$dprev/$scratchdir/53.0" "$d/$scratchdir"


    ### check if we completed at least one sp giving us a:
    #spcompleted_str='Max gradient component'
    #if [ $(grep -c "$spcompleted_str" "$d/$outfilename") -lt 1 ]; then

    #!# check if there is a HESS (BFGS guess) file (HESS should be newer than 132.0):
    if [ -f "$dprev/$scratchdir/HESS" ]; then
      #!# make hessian binary file:
      python ./write_fake_hessian.py "$dprev/$scratchdir/HESS"
      cp ./fake_hess_132.0 "$d/$scratchdir/132.0"
      rm ./fake_hess_132.0

    #!# check if there is a 132.0 file already there:
    elif [ -f "$dprev/$scratchdir/132.0" ]; then
      cp "$dprev/$scratchdir/132.0" "$d/$scratchdir"

    ## else we comment out read from hessian:
    else
      sed -i 's/GEOM_OPT_HESSIAN/!GEOM_OPT_HESSIAN/' "$d/$initinput"
    fi

    #!# make temp dir:
    cp -r "$d" ./tmp_${d}

    cd "$d"
    qsub -N "${count}_${dirprefix}" "$subscript" > job.id
    cd "$here"

  fi

  #!# else, calculation on queue or runnning:
  sleep "$sleeptime"

done







