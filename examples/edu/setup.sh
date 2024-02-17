#!/bin/bash

clear

[ "$SIZE" == "" ] && SIZE=4

COOK=cook$1slc

if [ "$NSLOTS" == "" ] ; then
  echo "serial version"
else
  COOK=${COOK}P1
  [ "$NSLOTS" -gt 2 ] && NSLOTS=2
  echo "parallel version (2 slots)"
fi

X=`ev "$SIZE*180+185"`
Y=`ev "$SIZE*120"`
export PLOTGEOMETRY="${X}x${Y}"
export RDFGGEOMETRY="${X}x${Y}"

X=`ev "$SIZE*150+185"`
Y=`ev "$SIZE*120"`
export SHOWGEOMETRY="${X}x${Y}"

X=`ev "$SIZE*180+185"`
Y=`ev "$SIZE*80"`
export LONGSHOW="${X}x${Y}"

X=`ev "$SIZE*180+185"`
Y=`ev "$SIZE*80"`
export SHOWCPGEOMETRY="${X}x${Y}"
