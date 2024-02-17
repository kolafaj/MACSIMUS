#!/bin/bash

source setup.sh ew
source slabname.sh

echo "%%%%%%%%% STEP A09 - preparing elongated box and melting a half of it %%%%%%%%%"
echo

if [ ! -e T.def ] ; then
  echo "ERROR! You likely have not entered temperature in Step 03. Repeat from Step 03."
  sleep 3
  exit
fi

if [ ! -e cryst.cfg ] ; then
  echo "ERROR! The expected files of Step A04 are missing. Repeat Step 04."
  sleep 3
  exit
fi

if [ ! -e box.def ] ; then
  echo "ERROR! File box.def with averaged box size from Step 04 is missing.  Repeat Step 04."
  sleep 3
  exit
fi

if [ -e $SLAB.cfg ] ; then
  echo "Step A09 has been already executed. For sure, I will erase everything and start again."
  echo
fi

if [ -e $SLAB.loc ] ; then
  echo "ERROR! The simulation of the same name is already running. I will try to interrupt it."
  echo "Wait a minute and try again."
  touch $SLAB.stp
  echo "[Enter]=continue"
  read
  exit
fi

cp cryst.cfg $SLAB.cfg

cat > $SLAB.def <<EOF
n=108*3               ! auxiliary variable (3x co predtim)
N[0]=n N[1]=n         ! the number of Na+ a Cl-

cutoff=8.607          ! elst cutoff [AA]
LJcutoff=cutoff       ! LJ cutoff   [AA]
!rdf.grid=20           ! radial distribution function histogram 20 bins/AA
el.epsk=2 el.epsr=0.4 ! electrostatic accuracy [K/AA]
el.diff=0.3           ! some warnings suppressedi
noint=30 h=0.1/noint  ! steps/cycle and step length [ps]
no=100                ! number of cycles
dt.plb=1              ! how often to record "playback" [ps]

thermostat="Andersen" ! random hits (Maxwell-Boltzmann distributon)

tau.T=1               ! thermostat time constant [ps] 

P=101325              ! pressure [Pa]
bulkmodulus=1e13/(T+300) ! bulk modulus estimate (for barostat) [Pa]
tau.P=1               ! barostat time constant [ps]
rescale="ZCM"         ! only the z-size of the box will be changed

init="start"          ! start from the previously recorded configuration, new m

! ERASE THIS AFTER THE FIRST STEP:
load.n[2]=3           ! replicate 3x in the direction of z
slab.T=2200           ! temperature in the central half of the box (mezi Tz0,Tz1)
z=(T-1300)/5000
slab.Tz0=0.25+z       ! range of z for the higher temperature slab.T (in L[2])
slab.Tz1=0.75-z       ! range of z for the higher temperature slab.T (in L[2])
T=800 tau.T=0.06      ! temperature thermostat time constant (faster)
tau.P=5               ! barostat time constant (faster)
EOF

echo "Now the box of Step A04 will be replicated in the z-direction and elongated."
echo "The central part of the box will be heated up to 1600 K, so that it will melt."
echo "The remaining part of the box will be cooled down to 1000 K."
echo "The x,y-sizes are not changed, the box may expand in the z-direction only;"
echo "the pressure in the z-direction is maintained by the barostat."

echo
echo "In Step A03, you have entered the following temperature:"
echo -n "   "
cat T.def
cat T.def >> $SLAB.def
echo "It will be used again."

echo
echo "In Step A04, the following box size was calculated at given T and p:"
echo -n "   "
cat box.def
cat box.def >> $SLAB.def 
BOX=`cut -c3-4 box.def`
if [ "$BOX" != 17 -a "$BOX" != 18 ] ; then
  echo "PROBLEM! The box size is not in the expected range."
  echo "[Enter]=continue, [Ctrl-C]=interrupt"
  read
fi

echo "This size will be used for the x,y-sizes of the box, the z-size will fluctuate."

cat >> $SLAB.def <<EOF
L[0]=x L[1]=y L[2]=3.2*x
z=T
T=800;
EOF

cat > $SLAB.get <<EOF
init="start"
dt.plb=0.5   ! show more often
tau.T=0.5 tau.rho=2 tau.P=0 no=10; ! short relaxation
L[2]=3.4*x; ! expand the box in the z-direction (do not melt yet)
L[2]=3.6*x; ! more expansion
slab.T=1600 ! temperature of the central slab (will melt)
T=1000      ! temperature of the other part (crystal kept)
tau.rho=0 tau.P=1  ! barostat
no=50; 
EOF

echo
echo "The calculations will be performed in the batch mode."
echo "[Enter]=continue, [Ctrl-C]=interrupt"
read

$COOK nacl $SLAB

echo
echo "Now I will show the slab melting, frames by 0.5 ps,"
echo "parallel projection and smaller spheres for clarity."
echo "[Enter]=continue, [Ctrl-C]=interrupt"
read

SHOWGEOMETRY=$LONGSHOW nice show $SLAB -s240 '-I%=rrrrrrrrrrrrri'$'\t'
