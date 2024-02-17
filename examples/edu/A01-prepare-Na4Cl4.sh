#!/bin/bash

source setup.sh ew

if [ -e na4cl4.* ] ; then
  echo "Erasing all files from previous Step A01..."
  sleep 1
  rm na4cl4.* {NA,CL}.* nacl.ble na4cl4.*
fi


cat > Na+.che <<EOF
Na+
parameter_set=sea
NAp1
EOF
cat > Cl-.che <<EOF
Cl+
parameter_set=sea
CLm1
EOF

# to show a cluster only
cat > Na4Cl4.che <<EOF
Na4Cl4
parameter_set=sea

NAp1 NAp1 NAp1 NAp1

CLn1 CLn1 CLn1 CLn1
EOF
blend -m0 Na4Cl4

if blend -o nacl NA CL ; then
  echo "-------------------------------------------------------------------"
  echo "I have prepared files with NaCl force field, atom sizes and colors."
  echo
  echo "FYI, the Lennard-Jones parameters follow:"
  echo "- Emin is the potential minimum, in kcal/mol"
  echo "- Rvdw = sigma*2^(1/6)/2, in AA"
  echo
  lemon Lennard-Jones +7 +h nacl.ble | lemon 'RvdW alpha[1-4]' '-  RvdW      alpha[1-4]  ' | cut -c1-80 | head -n12
  echo "[Enter]=continue, [Ctrl-C]=interrupt"
  read
else
  echo "ERROR! Force field generation failed, ask an expert."
  echo "Type [Enter]"
  read
  exit
fi

echo "%%%%%%%%%%% STEP A01 - preparing nanocrystal Na4Cl4 %%%%%%%%%%%"

echo "It's your turn!"
echo "The density of the model solid NaCl is 2.1 g/cm3, M(NaCl)=58.4 g/mol."
echo "Calculate the size L of the box containing Na4Cl4 in Angstrom (AA) and"
echo "enter it below (number only, not unit)."
echo "The value will be used to prepare atom coordinates of cluster Na4Cl4."
echo "Then, the cluster will be shown - click [quit] to continue."
echo

IF=2
while [ $IF -gt 0 ] ; do
  [ $IF == 1 ] && echo "This does not look good, try again!"
  echo -n "L = "
  read
  L=$REPLY
  L=${L/,/.}
  [ x$L == x ] && L=1
  IF=`ev "($L<5.5) | ($L>5.9)"`
done

echo "8 -3" > Na4Cl4.pla
echo "$L $L $L" >> Na4Cl4.pla
A=`ev "$L*.25"`
B=`ev "$L*.75"`
echo "$A $A $A" >> Na4Cl4.pla
echo "$B $B $A" >> Na4Cl4.pla
echo "$A $B $B" >> Na4Cl4.pla
echo "$B $A $B" >> Na4Cl4.pla
echo "$B $B $B" >> Na4Cl4.pla
echo "$A $A $B" >> Na4Cl4.pla
echo "$A $B $A" >> Na4Cl4.pla
echo "$B $A $A" >> Na4Cl4.pla

echo "The input file with Na4Cl4 coordinates is:"
cat Na4Cl4.pla

asc2plb Na4Cl4 -y

if [ ! -e Na4Cl4.plb ]; then
  echo "ERROR! The expected file is missing, run A01-prepare-Na4Cl4.sh again."
else
  echo "OK"
fi

show '-|O0.05' -Irr Na4Cl4 Na4Cl4.plb
