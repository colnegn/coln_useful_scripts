#!/bin/bash

## assume that we don't want to print headers:
print_headers_bool=0
while getopts ":h" opt; do
  case $opt in
    h)
      print_headers_bool=1
      shift
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


# Assume parameter passed in is a relative path to a directory.
# For brevity, we won't do argument type or length checking.

## loop over all arguments:
for a in $(seq 1 "$#"); do

  ## header:
  if [ $print_headers_bool -eq 1 ]; then
    echo "==> $1 <=="
  fi

  ## check if argument is file or directory:
  if [ -d $1 ]; then
    ## if directory, go to directory:
    ABS_PATH=`cd "$1"; pwd -P`              # double quotes for paths that contain spaces etc...
    echo "$ABS_PATH"
  else
    ## if not directory, go to containing directory:
    ABS_PATH=`cd "$(dirname $1)"; pwd -P`
    echo "$ABS_PATH/$(basename $1)"
  fi


  ## increment arguments:
  shift


  if [ $print_headers_bool -eq 1 ]; then
    ## only add extra newline if more arguments left:
    if [ $# -gt 0 ]; then
      echo
    fi
  fi

done


