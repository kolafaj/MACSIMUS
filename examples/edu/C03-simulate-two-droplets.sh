#!/bin/bash

source setup.sh ce
source Cdata.sh

echo "%%%%%%%%%%%%%% STEP C03 - two droplets from one + simulate %%%%%%%%%%%%%%"

cat > drop2.def <<EOF
n=$N
N[0]=n*2              ! the number of water molecules
rho=200               ! reference density [kg/m3]
cutoff=(n*6)^(1/3)    ! elst cutoff [AA]
x=12
cutoff=x+(cutoff-x)*(cutoff<x) ! max. 12
LJcutoff=10           ! Lennard-Jones cutoff
noint=5 h=0.01/noint  ! steps/cycle and step length [ps]
no=n                  ! number of cycles
dt.plb=1              ! how often to record "playback" [ps]
dt.prt=1              ! how often to print protocol [ps]
thermostat=0          ! no thermostat (adiabatic)
T=$T                 ! temperature [K]
tau.T=0.5             ! thermostat time constant [ps]
init="start"          ! start from the previously recorded configuration, new m
corr=0                ! no cutoff corrections
;
EOF

cat > drop2.get <<EOF
init="start"      ! use the previous configuration
load.n[0]=2       ! replicate it in the x-direction
no=0;             ! no simulation (just write)
EOF

cat drop2.get
cp drop.cfg drop2.cfg

if ! $COOK spce drop2 ; then
  echo ERROR
  exit
fi

cat > drop2.get <<EOF
init="start"
rho=200
x=(29915.272/rho*n)^(1/3)  ! box size
rho=0
L[0]=x*2                   ! box enlarged in x
L[1]=x*1.25                ! box enlarged in y
L[2]=x*1.25                ! box enlarged in z
load.L[0]=3                ! box scaled but not the coordinates (x)
load.L[1]=3                ! box scaled but not the coordinates (y)
load.L[2]=3                ! box scaled but not the coordinates (z)
shift[0]=1                 ! shift in x (to look better)
shift[1]=5                 ! shift in y
shift[2]=4                 ! shift in z
no=0;
EOF

if ! $COOK spce drop2 ; then
  echo "ERROR! The simulation failed."
  echo "Type [Enter]"
  read
  exit
fi

echo "--------------------------------------------------------------"
echo "The droplet has been replicated."
echo
echo "Both droplets will be assigned opposite velocities in a colliding course."
echo "The velocity should be given in m/s (number without unit)."
echo "The recommended value is 30, but you may try anything from 20 to 1000."
echo "In order to detect temperature increase on coelescence, use 25 to 35."
echo "Velocities greater than 2000 lead to vaporization."

echo -n "Enter the velocity: "
read
V=$REPLY

[ x$V = x ] && exit

cat > drop2.get <<EOF
init="start"
nshift=-n           ! half molecules right, half left
vshift[0]=$V/100    ! change of velocity [AA/ps]
dt.plb=1
rho=100             ! reference density
tau.T=0.2 thermostat="Berendsen"
no=10;
thermostat=0
init="append"
no=n*15;
EOF

echo ">>>>>> Simulation parameters:"
cat drop2.get

echo
echo "The simulation of droplet coalescence will be performed in the batch mode."
echo "It is possible to watch the running simulation (C04-watch-coalescence.sh)"
echo "(The stop request will be sent later.)"
echo "[Enter]=continue, [Ctrl-C]=interrupt"
read

NAME=`basename $PWD`
rm drop2.plb
$COOK spce drop2 &

sleep 3

echo
echo "Started in background, wait a minute and continue by Step C04"

echo "[Enter]=continue"
read
