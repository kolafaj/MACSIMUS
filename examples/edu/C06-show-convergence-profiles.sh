#!/bin/bash

source geometry.sh
echo "%%%%%%%%%%% STEP C05 - prubeh splyvani kapek %%%%%%%%%%%"
if [ ! -e drop2.cp ] ; then
  echo "ERROR - chybi soubor z kroku A04, simulaci nutno nejprve ukoncit."
  sleep 3
  exit
fi

rm *~

nice showcp -a -b10 -p4 drop2.cp Etot Tkin Epot

echo
echo "                       \\\\\\\\\\V/////"
echo "                      /  _    _  \\"
echo "                     |   O    O   |"
echo "                    @|     ..     |@"
echo "                      \\   \\__/   /"
echo "                       \\________/"
echo
