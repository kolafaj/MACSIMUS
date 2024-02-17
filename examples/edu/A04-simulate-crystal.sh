#!/bin/bash

source setup.sh ew

echo '%%%%%%%%%%% STEP A04 - productive simulation of crystal Na108Cl108 %%%%%%%%%%'

if [ ! -e cryst.cfg ] ; then
  echo "ERROR! File cryst.cfg is missing, repeat Step A03.sh"
  sleep 3
  exit
fi

if [ -e cryst.loc ] ; then
  echo "ERROR! The simulation of the same name is already running. I will try to interrupt it."
  int.sh cryst
  echo "Finished! Type [Enter] to continue."
  read
  exit
fi

echo ">>>>> simulation parameters:"
cat > cryst.get <<EOF
init="start"   ! start from the previous simulation
               ! measurements and recording initialized
dt.plb=1       ! how often to record a configuration (playback) [ps]
no=600;        ! the number of simulation cycles
EOF

cat cryst.get

echo
echo "Now the simulation is about to be started (batch mode)."
echo "[Enter]=continue, [Ctrl-C]=interrupt"
read

cat > cryst.cpi <<EOF
+Lx
EOF

[ -e cryst.mol ] && rm cryst.mol
$COOK nacl cryst

Lx=`staprt -nLx -ka -F%g cryst`
echo "x=$Lx y=$Lx" > box.def
echo
echo "The calculated averaged size of the simulation box (in AA):"
echo "  L=$Lx"
echo "It will be used later to set the box size."

echo "[Enter]=continue"
read
