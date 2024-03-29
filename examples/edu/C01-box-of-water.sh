#!/bin/bash

source setup.sh ce

echo $COOK

echo "The initial temperature will be 273 K"
echo "Enter the number of molecules in one droplet [200..400]:"
read
N=$REPLY

[ "$N" = "" ] && exit

cat > Cdata.sh <<EOF
# generated by C01-box-of-water.sh
T=273
N=$N
EOF

source Cdata.sh

echo "%%%%%%%%%%% STEP C01 - periodic box %%%%%%%%%%%"

blend -o spce -h SPCE || exit

cat > box.def <<EOF
n=$N
N[0]=n               ! the number of water molecules
rho=990              ! reference density [kg/m3]
cutoff=(n*3.7)^(1/3)  ! elst cutoff [AA]
x=12
cutoff=x+(cutoff-x)*(cutoff<x) ! max. 12
LJcutoff=cutoff       ! LJ cutoff   [AA]
noint=5 h=0.01/noint  ! steps/cycle and step length [ps]
no=n                  ! number of cycles = N
dt.plb=1              ! how often to record "playback" [ps]
thermostat="Berendsen"! thermostat type
T=$T                 ! temperature [K]
tau.T=0.5             ! thermostat time constant [ps]
corr=0                ! no cutoff corrections
init="start"          ! start from the previously recorded configuration
;
EOF

cat > box.get <<EOF
init="crystal"     ! crystal-like initial configuration
tau.T=0.15         ! thermostat time constant [ps]
T=$T-10
tau.P=0
dt.plb=0.1         ! record quite often [ps]
no=n/3;
EOF

echo "------- box.def --------"
cat box.def
echo "------- box.get --------"
cat box.get

echo
echo "Force field file spce.ble and global simulation parameters have been generated."
echo "Now the initial configuration will be constructed."
echo "Fast relaxation at density 990 kg/m3 follows."
echo "[Enter]=continue, [Ctrl-C]=interrupt"
read

rm box.cp
$COOK spce box

if [ ! -e box.cp ] ; then
  echo
  echo "ERROR! The simulation failed (missing box.cp)."
  echo "Type [Enter]"
  read
  exit
fi

showcp -p box P

