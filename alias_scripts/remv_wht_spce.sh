#!/bin/bash

if [ $# -ne 1 ]; then
  echo 'usage: 1.text string to remove bad \"name characters\"'
  exit
fi


#
#
#  MAKE A FUNCTION !!!!
#
#
#

inputstr=$1

wht_spce_str=$inputstr

echo " "

sansspace=`echo ${wht_spce_str} | sed -e 's/ /_/g'`
#echo $sansspace

sanscolon=`echo ${sansspace} | sed -e 's/:/_/g'`
#echo $sanscolon

sansequal=`echo ${sanscolon} | sed -e 's/=/_/g'`
#echo $sanscolon

sansqmark=`echo ${sansequal} | sed -e 's/?/_/g'`
#echo $sansqmark

sanscomma=`echo ${sansqmark} | sed -e 's/,/_/g'`
#echo $sanscomma

sansfslash=`echo ${sanscomma} | sed -e 's#/#__#g'`
#echo $sansfslash

sanshash=`echo ${sansfslash} | sed -e 's/#/_/g'`
#echo $sanshash

sansperiod=`echo ${sanshash} | sed -e 's/\./_/g'`
#echo $sansperiod

sansplus=`echo ${sansperiod} | sed -e 's/+/-POS__/g'`
#echo $sansplus

sanslparen=`echo ${sansplus} | sed -e 's/(/--/g'`
#echo $sanslparen

sansrparen=`echo ${sanslparen} | sed -e 's/)/--/g'`
#echo $sansrparen

sansarrow=`echo ${sansrparen} | sed -e 's/→/-to-/g'`
#echo $sansarrow

sansdagger=`echo ${sansarrow} | sed -e 's/†/_/g'`
#echo $sansdagger

sansdot=`echo ${sansdagger} | sed -e 's/·/-/g'`
#echo $sansdot

sansasterisk=`echo ${sansdot} | sed -e 's/*/-/g'`
#echo $sansasterisk

sansdash=`echo ${sansasterisk} | sed -e 's/—/-/g;s/—/-/g'`
echo $sansdash

#sansbslash=`echo ${sansperiod} | sed -e 's#\\#_#g'`
#echo $sansbslash

echo " "

