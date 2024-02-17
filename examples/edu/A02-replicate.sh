#!/bin/bash

source setup.sh ew
echo "%%%%%%%%%%% STEP A02 - replicate Na4Cl4 -> Na108Cl108 %%%%%%%%%%%"

if [ ! -e na4cl4.plb ] ; then
  echo "ERROR! File na4cl4.plb is missing, repeat Step A01"
  sleep 3
  exit
fi

if [ -e cryst300.def ] ; then
  echo "Erasing all files from previous Step A02..."
  sleep 1
  rm cryst300.*
fi

echo
echo "I will replicate the crystal from Step A01 3x in each direction."
echo "1 ps long simulation with the Andersen thermostat will follow."
echo "[Enter]=continue, [Ctrl-C]=interrupt"
read
cat > cryst300.def <<EOF
n=108                 ! auxiliary variable
N[0]=n N[1]=n         ! the number of Na+ a Cl-

rho=2050              ! reference density [kg/m3]

cutoff=8.607          ! elst cutoff [AA]
LJcutoff=cutoff       ! LJ cutoff   [AA]
rdf.grid=20           ! radial distribution function histogram 20 bins/AA
el.epsk=2 el.epsr=0.4 ! electrostatic accuracy [K/AA]
el.diff=0.3           ! some warnings suppressedi
el.Perr=1e7           ! some warnings suppressedi
noint=30 h=0.1/noint  ! steps/cycle and step length [ps]
no=100                ! number of cycles
dt.plb=1              ! how often to record "playback" [ps]

thermostat="Andersen" ! random hits (Maxwell-Boltzmann distributon)
T=300                 ! temperature [K]
tau.T=1               ! thermostat time constant [ps] 

P=101325              ! pressure [Pa]
bulkmodulus=2e13/(T+300) ! bulk modulus estimate (for barostat) [Pa]
tau.P=2               ! barostat time constant [ps]

init="start"          ! start from the previously recorded configuration, new m

! ERASE THIS AFTER THE FIRST STEP:
load.n[0]=3           ! replicate 3x in the direction of x
load.n[1]=3           ! replicate 3x in the direction of y
load.n[2]=3           ! replicate 3x in the direction of z
;
EOF
cat > cryst300.get <<EOF
init="plb" tau.P=0 tau.rho=1
dt.plb=0.1 no=10;
EOF

echo

if ! $COOK nacl cryst300 na4cl4.plb ; then
  echo "ERROR! The simulation has crashed."                       
  rm cryst300.cfg
fi

echo
if [ ! -e cryst300.cfg ] ; then
  echo "ERROR! The expected output is missing, repeat Step A02"
else
  echo "The simulation has stopped, now I will show the trajectory."
  echo "[Enter]=continue, [Ctrl-C]=interrupt"
  read

  nice show cryst300 '-I%i---'
fi
