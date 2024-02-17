#!/bin/bash

source setup.sh ce

echo "%%%%%%%%%%% STEP B01 - select the solute, initialize, NVT relaxation %%%%%%%%%%%"

N=200
T=300
C=9
nfactor=1

for (( ;; )) ; do
  echo "Enter molecule name from the list."
  echo "Noble gases (Lennard-Jones):"
  echo "  He"
  echo "  Ne"
  echo "  Ar"
  echo "  Kr"
  echo "  Xe"
  echo "Ions: (Lennard-Jones + charge):"
  echo "  Li+"
  echo "  Na+"
  echo "  Cl-"
  echo "  Ca++"
  echo "Endofullerenes (atom in C60, SLOW):"
  echo "  Xe@C60"
  echo "  Ca++@C60"
  echo "  Cl-@C60"
  read
  M=$REPLY
  [ -e $M.che ] && break

  echo "ERROR! \'$M.che\' does not exist - try again."
  sleep 1
done

if [[ $M == *"@"* ]] ; then 
  N=350
  C=11
  nfactor=2
fi

cp $M.che solute.che

rm $M.?ol $M.plb

if ! blend -o ff $M -h SPCE ; then 
  echo "ERROR! No force field has been generated."
  echo "type [Enter]"
  read
  exit
fi

cat > solution.def <<EOF
N[0]=1                ! 1 solute
N[1]=$N              ! the number of water molecules SPC/E
rho=1000              ! reference density [kg/m3]
cutoff=$C              ! elst cutoff [AA]
LJcutoff=cutoff       ! LJ cutoff   [AA]
rdf.grid=10           ! bins per 1 AA for the radial distribution function histogram
noint=6 h=0.01/noint  ! steps/cycle and step length [ps]
no=300                ! number of cycles
dt.plb=1              ! how often to record "playback" [ps]
thermostat="Berendsen"! thermostat type
T=$T                 ! temperature [K]
tau.T=0.5             ! thermostat time constant [ps]
P=101325              ! pressure [Pa]
bulkmodulus=2e9       ! bulk modulus estimate (for barostat) [Pa]
tau.P=10              ! barostat time constant
center.cmn=1          ! central force for the 1st molecule
center.cmK[0]=300     ! x-force in K/AA
center.cmK[1]=300     ! y-force in K/AA
center.cmK[2]=300     ! z-force in K/AA
pins=0                ! do not adjust box size with init="random"
el.epsq=9e9           ! will accept charged system
a=$nfactor            ! used to double 'no' for endofullerenes
;
EOF

echo "------ solution.def ------"
cat solution.def

cat > solution.get <<EOF
init="random" ! random initial configuration
tau.T=0.15    ! fast cooling
T=$T-10      ! bit lower temperature
tau.P=0       ! barostat off
no=150*a;     ! the number of cycles by 0.01 ps
EOF

echo "------ solution.get ------"
cat solution.get

echo
echo "Force field file ff.ble and global simulation parameters have been generated."
echo "Now the initial configuration will be constructed."
echo "The solute is attracted to the box center."
echo "Fast relaxation at density 1000 kg/m3 follows."
echo "Type [Enter]"

read

rm solution.cp solution.mol

$COOK ff solution

if [ ! -e solution.cp ] ; then
  echo
  echo "ERROR! The simulation failed (missing solution.cp)."
  echo "Type [Enter]"
  read
  exit
fi

echo ">>>> convergence profiles <<<<<<"

nice showcp -p5 solution -2 -3 -5
