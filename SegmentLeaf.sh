#!/bin/bash

do_run()
{
  echo $*
  if $* ; then
    true
  else
    exit 1
  fi
}




checkpython()
{
  if [ -z "$PYTHONPATH" ]
    then
    export PYTHONPATH=/usr/bin/python
  fi
}



callSegment()
{
  #Focus=$1
  #Nums=$2
  #Pics=$3
  echo $Pics
exit 0
  do_run python segTestAra.py -f ${Focus} -n ${Nums} -p ${Pics}
}


Focus=''
Nums=''
Pics='hola'

while getopts :f:n:p:o opt; do
   case "$opt" in
     f) Focus="$OPTARG";;	
     n) Nums="$OPTARG";;	
     f) Pics="$OPTARG";;	
     o) isdef=1;;	
     \?) help_ani;;
   esac
done	

echo $Pics
checkpython
callSegment

