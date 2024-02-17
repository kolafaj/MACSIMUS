#!/bin/bash

source setup.sh ew

echo "%%%%%%%%%%% STEP A08 - watching the results of Step A07 %%%%%%%%%%%"
if [ ! -e box.def ] ; then
  echo "ERROR! The expected files of Step A07 are missing. Repeat Step 07."
  sleep 3
  exit
fi

rm *~

./res.sh melt rho
