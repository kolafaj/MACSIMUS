#!/bin/bash

source setup.sh ew

echo "%%%%%%%%% STEP A03 - equilibrating crystal Na108Cl108 at given T and p %%%%%%%%"

if [ ! -e cryst300.cfg ] ; then
  echo
  echo "ERROR! File cryst300.cfg is missing, repeat Step A02"
  sleep 3
  exit
fi

if [ -e cryst.cfg ] ; then
  echo "Step A03 has been already executed."
  if ./yes.sh "Use the configuration from Step A02 and repeat all the initialization (y/N)?" ; then
    cp cryst300.cfg cryst.cfg
    cat > cryst.get <<EOF
init="start"
dt.plb=1
no=100;
EOF
  fi
else
  cp cryst300.cfg cryst.cfg
  cat > cryst.get <<EOF
init="start"
dt.plb=1
no=100;
EOF
fi
  
while [ x$R == x ] ; do

  echo
  echo "Everything is ready to simulate a crystal."

  while [ "$T" == "" ] ; do
    echo -n "Enter the temperature (in K) in the interval about 1200-1400: "
    read
    T=$REPLY
  done

  T=${T/K/}
  T=${T// /}

  echo "T=$T" > T.def
  echo "SLAB=slab_${T}K" > slabname.sh

  cat > cryst.def <<EOF
n=108                 ! auxiliary variable
N[0]=n N[1]=n         ! the number of Na+ a Cl-

rho=1900              ! reference density [kg/m3]

cutoff=8.607          ! elst cutoff (for the Ewald summation)[AA]
LJcutoff=cutoff       ! Lennard-Jones cutoff [AA]
rdf.grid=20           ! bins per 1 AA for the radial distribution function histogram
el.epsk=2 el.epsr=0.4 ! electrostatic accuracy [K/AA]
el.diff=0.3           ! some warnings suppressed
noint=30 h=0.1/noint  ! steps/cycle and step length [ps]
no=100                ! number of cycles
dt.plb=1              ! how often to record "playback" [ps]

thermostat="Andersen" ! random hits (Maxwell-Boltzmann distributon)
T=$T                 ! temperature [K]
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
  cat cryst.def
  echo ">>>>> more simulation parameters:"
  cat cryst.get
  echo "[Enter]=continue, [Ctrl-C]=interrupt"
  read

  $COOK nacl cryst

  if [ ! -e cryst.cp ] ; then
    echo "ERROR! The expected output is missing, repeat Step A03"
    echo "(spatne zadana temperature?)"
    echo "Type [Enter]..."
    read
    exit
  fi

  echo
  echo "==============================================================================="
  echo "Now I will show the time dependence of temperature, energy, and density."
  echo "Hint: you may close all graphs at once by button [kill all]"
  echo "If the graphs do not converge, repeat equilibration!"
  echo "[Enter]=continue, [Ctrl-C]=interrupt"
  read
  nice showcp -p5 cryst Tkin Epot rho

  cat > cryst.get <<EOF
init="append"
dt.plb=1
no=100;
EOF

  echo
  if ./yes.sh "Repeat equilibrating with init=\"append\" (Y/n)?" y ; then R="" ; else R=n ; fi

done
