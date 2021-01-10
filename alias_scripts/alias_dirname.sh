#!/bin/bash


##################################################
###  getopts  ####################################


## defaults:
num_dirs=1
full_path_bool=0


## number of input arguments to shift after getopts:
num_shift=0
while getopts ":n:f:" opt; do
  case $opt in
    n)
      num_dirs=$OPTARG
      num_shift=$((num_shift+2))
      ;;
    f)
      full_path_bool=1
      num_shift=$((num_shift+2))
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done

## now shift input arguments:
shift $num_shift

#================================================#
#================================================#



##################################################
###  inputs & usage  #############################

if [ $# -ne 1 ]; then
  echo '
  
         usage:

  options:    
              -n <number of (containing) directories to look back toward>
                 (
                  e.g.: 
                         ...]$ dname -n 3 ./dir1/dir2/dir3/dir4/dir5/dir6/file1
                         >>> dir4

                   e.g.:
                         ...]$ dname -n 3 ./dir1/dir2/dir3/dir4/dir5/dir6
                         >>> dir3
                 )

              -f <full path bool>
                 (
                  i.e.:
                         0 ~ just directory name
                         1 ~ full path
                 )

  arguments:
    1.path to extract directory name from

          '
  exit
fi

here=`pwd`

in_path=$1

#================================================#
#================================================#



##################################################
###  call dirname function  ######################

function dname {

  local in_dpath=$1
  local count=$2

  ## found the desired directory:
  if [ $count -eq 0 ]; then

    ## return "full" path:
    if [ $full_path_bool -eq 1 ]; then
      printf "$in_dpath"
    else
      printf "$(basename $in_dpath)"
    fi

  ## have not yet found directory, so recursively-call self:
  else
  
    ## "increment" counter:
    count=$((count - 1))

    ## call self:
    out_dname=$(dname $(dirname $in_dpath) $count)


    ## return value:
    printf $out_dname

  fi

  }


#echo $out_dname

#================================================#
#================================================#




#!# call function:
echo "$(dname $in_path $num_dirs)"




exit

