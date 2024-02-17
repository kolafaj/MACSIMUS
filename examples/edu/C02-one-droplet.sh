#!/bin/bash

source setup.sh ce

source Cdata.sh

echo "%%%%%%%%%%% STEP C02 - a droplet by enlarging the box  %%%%%%%%%%%"

cp box.cfg drop.cfg

cat > drop.def <<EOF
n=$N
N[0]=n               ! the number of water molecules
rho=990               ! reference density [kg/m3]
cutoff=(n*3.7)^(1/3)  ! elst cutoff [AA]
x=12
cutoff=x+(cutoff-x)*(cutoff<x) ! max. 12
LJcutoff=cutoff       ! LJ cutoff   [AA]
noint=5 h=0.01/noint  ! steps/cycle and step length [ps]
no=300                ! number of cycles
dt.plb=1              ! how often to record "playback" [ps]
thermostat="Berendsen"! thermostat type
T=$T                 ! temperature [K]
tau.T=0.5             ! thermostat time constant [ps]
corr=0                ! no cutoff corrections
init="start"          ! start from the previously recorded configuration
;
EOF

  echo ">>>>> global simulation parameters:"
cat > drop.get <<EOF
init="start"
rho=200             ! density [kg/m3]
cutoff=(n*6)^(1/3)  ! elst cutoff [AA]
x=12                ! auxiliary variable
cutoff=x+(cutoff-x)*(cutoff<x) ! cutoff max. 12
LJcutoff=cutoff
load.L[0]=3         ! use the calculated box x-size
load.L[1]=3         ! use the calculated box y-size
load.L[2]=3         ! use the calculated box z-size
a=cutoff*2/3
shift[0]=a          ! shift the whole configuration in x
shift[1]=a          ! shift the whole configuration in y
shift[2]=a          ! shift the whole configuration in z
tau.T=0.3
no=n*2;
EOF

  echo ">>>>> global simulation parameters:"
  cat drop.def
  echo ">>>>> morel simulation parameters:"
  cat drop.get

echo
echo "The box will be enlarged, the molecules shifted"
echo "[Enter]=continue, [Ctrl-C]=interrupt"
read
$COOK spce drop

if [ ! -e drop.cp ] ; then
  echo
  echo "ERROR! The simulation failed (missing drop.cp)."
  echo "Type [Enter]"
  read
  exit
fi

echo ">>>> look at the configuration and the convergence profiles <<<<<<"
echo "NB: RDF has not been calculated."

res.sh drop
