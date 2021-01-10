#!/bin/bash


#!1!# hard-coded sleep time for while-loop:
sleep_time='5h'
sleep_time_int=5

#!2!# HARD-CODED maximum number of loop iterations:
max_loop=48



if [ $# -ne 3 ]; then
  echo "usage: 1.process id number for job to monitor; 2.recognizable name for job; 3.email address to contact when job is done"
  exit
fi

pid=$1
job_name=$2
email_address=$3


## send out welcome email:
echo "
You are now subscribed to email_when_done.sh notifications for job ${job_name}.


Expect another email within $((sleep_time_int * max_loop)) hours, regardless of whether the job has finished running.


" | mail -s "Your ${job_name} job is running" $email_address



## check status of job with ps:
#job_status=$(ps -u $USER | grep -c "$pid")
job_status=$(ps -A | grep -c "$pid")


#!2!# don't want to loop forever:
loop_counter=0


## loop while job-id is in processes:
while [ $job_status -ne 0 ]; do

  #!1!# sleep:
  sleep $sleep_time

  ## check status of job for next iteration:
  job_status=$(ps -u $USER | grep -c "$pid")


  #!2!# check if the while-loop reached maximum iterations:
  loop_counter=$((loop_counter + 1))
  if [ $loop_counter -ge $max_loop ]; then
    email_str="$(date)    $job_name (pid: ${pid}) is still running after $max_loop x $sleep_time"
    echo "$email_str" | mail -s "$job_name (pid: ${pid}) is STILL RUNNING" $email_address
    exit
  fi

done


#!# if we made it here, job must be finished:
email_str="$(date)    $job_name (pid: ${pid}) finished within the past $sleep_time"
echo "$email_str" | mail -s "$job_name (pid: ${pid}) is COMPLETE" $email_address

exit



