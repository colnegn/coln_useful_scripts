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

































































trap '' 2  
trap '' 20

nits=40
space='  '
sleeptime='1.5'
invacceleration='0.92'
swoosh_timer=$((nits / 2))
for i in `seq 1 $nits`; do
  
  sleeptime=$(echo "$sleeptime * ($invacceleration)" | bc -l)

  sleep ${sleeptime}s
  clear

  echo "

                                            !!!!!!!!!!!!!!!!!!!!!!!!
                                            ! RESISTANCE IS FUTILE !
                                            !!!!!!!!!!!!!!!!!!!!!!!!

  "
  


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

echo "

                                            !!!!!!!!!!!!!!!!!!!!!!!!
                                            ! RESISTANCE IS FUTILE !
                                            !!!!!!!!!!!!!!!!!!!!!!!!

"

echo "$bird"

sleep 3s
clear

echo "

                                            !!!!!!!!!!!!!!!!!!!!!!!!
                                            ! RESISTANCE IS FUTILE !
                                            !!!!!!!!!!!!!!!!!!!!!!!!

"

echo "$caw"
sleep 2s


trap 2
trap 20






















































