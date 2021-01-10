#!/bin/bash


############################################################
###  defaults  #############################################

## default comment character:
com_char='#'
## default number of characters in header-comment:
len_header=50
## default number of left-padding characters:
num_left_pad=3
## default bottom-header character:
bot_head_char='='
## default number of 'extra' lines (for big headers, etc.):
extra_lines=0
## array for extra line characters:
declare -a extra_char
## default 'suppress bottom-header lines bool:
suppress_bot_head=0

#==========================================================#
#==========================================================#


############################################################
###  usage  ################################################

if [ $# -eq 0 ]; then
  echo "
  
        usage: 
  
  options: 
           -c  < comment character > (needs escape char!)
           -l  < length of header >
           -p  < number of padding characters on left >
           -b  < bottom \'header\' character > (needs escape char!)
           -s  suppress bottom \'header\'...

  arguments:
    input a string (spaces are fine) to get a nice header-comment
 

  default comment char:                   $com_char
  default number of characters in header: $len_header
  
  "
  exit
fi

#==========================================================#
#==========================================================#



############################################################
###  getopts  ##############################################

## number of input arguments to shift after getopts:
num_shift=0

while getopts ":c:l:p:b:e:s" opt; do
  case $opt in
    c)
      com_char=$OPTARG
      num_shift=$((num_shift+2))
      ;;
    l)
      len_header=$OPTARG
      num_shift=$((num_shift+2))
      ;;
    p)
      num_left_pad=$OPTARG
      num_shift=$((num_shift+2))
      ;;
    b)
      bot_head_char=$OPTARG
      num_shift=$((num_shift+2))
      ;;
    e)
      ## characters for extra lines (0 for header, 1 for bot-header)
      extra_char[$extra_lines]=$OPTARG
      ## increment number of extra lines (1 for header, 2 for both)
      extra_lines=$((extra_lines + 1))
      num_shift=$((num_shift+2))
      ;;
    s)
      ## 1 ~ suppress; 0 ~ don't suppress
      suppress_bot_head=1
      num_shift=$((num_shift+1))
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

#==========================================================#
#==========================================================#



############################################################
###  get header from remaining input arguments  ############

full_str="$1"
shift

for i in `seq 1 $#`; do
  full_str="${full_str} $1"
  shift
done

## number of characters in full input string:
len_str=${#full_str}

## want to exit if input string is longer than $((len_header - 7)):
if [ $((len_str + 6)) -ge $len_header ]; then
  echo "ERROR: input string is too long... exiting..."
  exit
fi

#==========================================================#
#==========================================================#



############################################################
###  generate header  ######################################

## first line of header-comment is only the comment character:
head1_line=''
for c in `seq 1 $len_header`; do
  head1_line="${head1_line}${com_char}"
done


## second line of header-comment has string with comments on either side:
head2_line=''
for c in `seq 1 $num_left_pad`; do
  head2_line="${head2_line}${com_char}"
done

## now insert (input) string into $head2_line with 2 spaces on either side:
head2_line="${head2_line}  ${full_str}  "

## now finish $head2_line with remaining number of comment characters:
for c in `seq ${#head2_line} $((len_header - 1))`; do
  head2_line="${head2_line}${com_char}"
done


## extra comment line for header:
if [ $extra_lines -ge 1 ]; then
  extra_header="${com_char}"
  for c in `seq 1 $((len_header - 2))`; do
    extra_header="${extra_header}${extra_char[0]}"
  done
  extra_header="${extra_header}${com_char}"

  ## print header (with extra line) to stdout:
  echo ""
  echo "$head1_line"
  echo "$extra_header"
  echo "$head2_line"
  echo ""

## else no extra comment line:
else
  ## print header to stdout:
  echo ""
  echo "$head1_line"
  echo "$head2_line"
  echo ""
fi

#==========================================================#
#==========================================================#



############################################################
###  generate bottom-header  ###############################

if [ $suppress_bot_head -eq 0 ]; then

  bot_line="${com_char}"
  for c in `seq 1 $((len_header-2))`; do
    bot_line="${bot_line}${bot_head_char}"
  done
  bot_line="${bot_line}${com_char}"


  ## extra comment line for bottom-header:
  if [ $extra_lines -eq 2 ]; then
    extra_bot="${com_char}"
    for c in `seq 1 $((len_header - 2))`; do
      extra_bot="${extra_bot}${extra_char[1]}"
    done
    extra_bot="${extra_bot}${com_char}"
  
    ## print header (with extra line) to stdout:
    echo ""
    echo "$bot_line"
    echo "$extra_bot"
    echo "$bot_line"
    echo ""
  
  ## else no extra comment line:
  else
    ## print bottom-header to stdout:
    echo ""
    echo "$bot_line"
    echo "$bot_line"
    echo ""
  fi
fi 

#==========================================================#
#==========================================================#


exit

