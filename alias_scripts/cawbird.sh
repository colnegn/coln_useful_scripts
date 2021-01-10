#!/bin/bash





















































caw='















              ________________
             /                \
             |  CAW !!!!!!!   |
             |                |
              \  ____________/
        __     \/  __    
       /  \  oCo  /  \
      /    \[___]/    \
     /     / | | \     \
     \/\/\/  w w  \/\/\/



'





bird='                
















                          
                          
                          
                          
        __         __     
       /  \  oCo  /  \    
      /    \[___]/    \   
     /     / | | \     \  
     \/\/\/  w w  \/\/\/  



'


flapbird='                
                          
                          

















                          
                          
     /\/\/\       /\/\/\  
     \     \ oCo /     /  
      \____/[___]\____/   
             | |          
             w w          


'




swooshbird='                
                          
                          












               
   ===========
   SWOOSH MODE
   ===========
                                       
                          
       ___   \ /   __     
----- /   \  oCo  /  \    
 --- /     \[___]/    \   
  --/     /  / / \     \  
    \/\/\/  w w   \/\/\/  



'













































SCRIPTPATH="$( cd "$(dirname "$0")" ; pwd -P )"
trap_str="bash $SCRIPTPATH/noresistcawbird.sh"

trap "$trap_str" 2
trap "$trap_str" 20

#trap 'bash ~/alias_scripts/noresistcawbird.sh' 2  
#trap 'bash ~/alias_scripts/noresistcawbird.sh' 20

nits=40
space='  '
sleeptime='1'
invacceleration='0.9'
swoosh_timer=$((nits / 2))
for i in `seq 1 $nits`; do
  
  sleeptime=$(echo "$sleeptime * ($invacceleration)" | bc -l)

  sleep ${sleeptime}s
  clear

  bird="$(echo "$bird" | sed "s/^/${space}/g")"
  flapbird="$(echo "$flapbird" | sed "s/^/${space}/g")"
  swooshbird="$(echo "$swooshbird" | sed "s/^/${space}/g")"
  caw="$(echo "$caw" | sed "s/^/${space}/g")"

  if [ $i -le $swoosh_timer ]; then
    if [ $(echo "$i % 2" | bc) -eq 1 ]; then
      echo "$bird"
    else
      echo "$flapbird"
    fi
  else
    echo "$swooshbird"
  fi

done

clear
echo "$bird"

sleep 3s
clear
echo "$caw"



trap 2
trap 20




























































