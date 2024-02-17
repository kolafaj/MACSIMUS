#!/bin/bash

source setup.sh ce

echo "%%%%%%%%%%% STEP B03 - simulation in the NPT ensemble %%%%%%%%%%%"

if [ ! -e solution.cfg ] ; then
  echo
  echo "ERROR! File solution.cfg is missing, repeat Step B01.sh"
  sleep 3
  exit
fi

cat > solution.get <<EOF
init="start"
tau.P=5
no=1000;
EOF

echo ">>>>>> simulation parameters:"
cat solution.get

rm solution.cp
echo
echo "Productive run of a solute in water."
echo "Type [Enter] to simulate, Ctrl-C to interrupt."
read

echo $COOK ff solution                                                                                              

$COOK ff solution

if [ ! -e solution.cp ] ; then
  echo "ERROR! File solution.cp is missing."
  echo "Type [Enter]..."
  read
  exit
fi
