#!/bin/bash

source setup.sh ew

echo "%%%%%%%%%%% STEP A06 - crystal melting %%%%%%%%%%%"

echo
if [ -e melt.cfg ] ; then

  echo "Step A06 has been already executed."
  if ./yes.sh "Use the configuration from Step A02 and repeat the initialization (y/N)?" ; then
    cp cryst300.cfg melt.cfg
    cat > melt.get <<EOF
init="start"
dt.plb=1
no=200;
EOF
  fi
else
  cp cryst300.cfg melt.cfg
  cat > melt.get <<EOF
init="start"
dt.plb=1
no=200;
EOF
fi

while [ x$R == x ] ; do

  if [ ! -e cryst300.cfg ] ; then
    echo "ERROR! File cryst300.cfg is missing, repeat Step A02"
    sleep 3
    exit
  fi

  echo
  echo "Steps A06-A08 (melting) can be skipped."
  echo "To do so, type Ctrl-C and continue by Step A09"
  echo

  T=1900

 while [  "$T" == "" ] ; do
   echo "Select temperature (in K) in interval 1800-2500"
   read
   T=$REPLY
 done

  cat > melt.def <<EOF
n=108                 ! auxiliary variable
N[0]=n N[1]=n         ! the number of Na+ a Cl-

rho=1700              ! reference density [kg/m3]

cutoff=8.607          ! elst cutoff [AA]
LJcutoff=cutoff       ! LJ cutoff   [AA]
rdf.grid=20           ! radial distribution function histogram 20 bins/AA
el.epsk=2 el.epsr=0.4 ! electrostatic accuracy [K/AA]
el.diff=0.3           ! some warnings suppressed
el.Perr=1e7           ! some warnings suppressed
noint=30 h=0.1/noint  ! steps/cycle and step length [ps]
no=100                ! number of cycles
dt.plb=1              ! how often to record "playback" [ps]

thermostat="Andersen" ! random hits (Maxwell-Boltzmann distributon)
T=$T                ! temperature [K]
tau.T=1               ! thermostat time constant [ps]

P=101325              ! pressure [Pa]
bulkmodulus=1e13/(T+300) ! bulk modulus estimate (for barostat) [Pa]
tau.P=2               ! barostat time constant [ps]

init="start"          ! start from the previously recorded configuration,
                      ! measurement and recording initialized
;
EOF

  echo
  echo ">>>>> global simulation parameters:"
  cat melt.def
  echo ">>>>> more simulation parameters:"
  cat melt.get
  echo "The simulation will be performed on the server."
  echo "[Enter]=simulate, [Ctrl-C]=interrupt"
  read

  $COOK nacl melt

  if [ ! -e melt.cp ] ; then
    echo "ERROR! The expected output is missing, repeat Step A03"
    exit
  fi

  echo
  echo "==============================================================================="
  echo "Now I will show the time dependence of temperature, energy, and density."
  echo "[Enter]=show, [Ctrl-C]=interrupt"
  read
  nice showcp -p5 melt Tkin Epot rho

  cat > melt.get <<EOF
init="append"
dt.plb=1
no=300;
EOF

  if ./yes.sh "Repeat equilibrating with init=\"append\" (Y/n)?" y ; then R="" ; else R=n ; fi

done
