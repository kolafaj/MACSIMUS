#!/bin/bash

source setup.sh ew
source slabname.sh

echo "%%%%%%%%%%% STEP A11 - progress of melting %%%%%%%%%%%"

if [ ! -s $SLAB.plb ] ; then
  echo "ERROR! File $SLAB.plb not (yet) available. Wait and re-try."
  echo "[Enter]=continue"
  read
  exit
fi

SHOWGEOMETRY=$LONGSHOW nice show $SLAB -s250 '-I%=rrrrrrrrrrrrri'$'\t'

if [ -e $SLAB.loc ] ; then
  if ./yes.sh "The slab simulation is still running. Interrupt (Y/n)?" y ; then
    int.sh $SLAB
  fi
fi

echo "Removing backup files..."
rm *~
sleep 1
