#!/bin/bash


#!# default values:
Nrows_down=0
Ncolumn=0
use_awk_bool=1
delimiter=' '
Noccurrence=1



num_shift=0

while getopts ":r:c:d:n:" opt; do
  case $opt in
    r)
      Nrows_down="$OPTARG"
      num_shift=$((num_shift+2))
      ;;
    c)
      Ncolumn="$OPTARG"
      num_shift=$((num_shift+2))
      ;;
    d)
      delimiter="$OPTARG"
      ## if we have a delimiter, don't use awk:
      use_awk_bool=0
      num_shift=$((num_shift+2))
      ;;
    n)
      Noccurrence="$OPTARG"
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


## usage:
if [ $# -ne 2 ]; then
  echo "
usage:
  argument:
    1.grep string
    2.file name to grep within

  options:
    -r number of rows to look down (argument of grep -A; use '0' for same line as grep string)
    -c column number (for awk or cut)
    -d delimiter character (turns off awk, turns on cut)
    -n number of occurrence to grab
"
  exit
fi

## grep string from required argument:
grep_str="$1"
shift
## file name from required argument:
fname="$1"
shift




#!# determine line number with Noccurrence:
bot_line="$((Noccurrence * (Nrows_down + 2) - 1))"
#                                      ^^^^^^^^ grep -A adds a '--' line between each result

#!# do it:
if [ $use_awk_bool -eq 1 ]; then
  ## with awk:
  grep -A"$Nrows_down" "$grep_str" "$fname" | head -$bot_line | tail -1 | awk -v n="$Ncolumn" '{print $n}'
else
  ## with cut:
  grep -A"$Nrows_down" "$grep_str" "$fname" | head -$bot_line | tail -1 | cut -d"$delimiter" -f"$Ncolumn"
fi


exit


