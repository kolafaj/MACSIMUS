#!/bin/bash

source setup.sh ce

echo "%%%%%%%%%%% STEP B05 - solvation shell %%%%%%%%%%%"

if [ ! -e solution.cfg ] ; then
  echo
  echo "ERROR! File solution.cfg is missing, repeat Step B01"
  sleep 3
  exit
fi

cat <<EOF
The solvation shell of radius 5 AA around the central atom will be shown.
The shell thickness can be changed by hot keys jJ (there is no mouse button).
Showing hydrogen bonds is turned on by button [H-bonds] or Ctrl-H.
The H-bond length can be changed by buttons [+] [-].
Note that the length corresponds to the first minimum on th H-O RDF.
Watch the orientation of water molecules near the solute!
EOF

nice show -l -j5 -Jb -I@ solution
