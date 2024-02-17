#!/bin/bash
if [ $# -le  0 ]; then
  echo "Call by:"
  echo "  nmf.sh MOLNAME[.che] [BLEND-OPTIONS]"
  echo "Examples:"
  echo "  nmf.sh biphenyl -M.4"
  echo "  nmf.sh biphenyl.che -r2"
  echo "See also:"
  echo "  Aharm"
  exit
fi

if [ "${BLENDPATH}" == "" ] ; then
  export BLENDPATH=/home/$USER/macsimus/blend/data
fi
N=${1/%.che/}
if [ -e $N.plb -o -e $N.mol ] ; then
  echo "$N.plb or $N.mol exist. Erase before run (y/N)?"
  read
 if [ $REPLY = "y" ] ; then
   rm $N.plb $N.mol 
 fi
fi

blend -N -M-1 -P999 -m99999 -p -g10 $2 $1
shownmf.sh $N.nmf
