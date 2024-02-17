#!/bin/bash

source setup.sh ew
source slabname.sh

echo "%%%%%%%%%%% STEP A10 - slab simulation %%%%%%%%%%%"

if [ ! -e $SLAB.cfg ] ; then
  echo
  echo "ERROR! File $SLAB.cfg is missing, repeat Step A09"
  sleep 3
  exit
fi

if [ -e $SLAB.loc ] ; then
  echo "ERROR! The simulation of the same name is already running. I will try to interrupt it."
  int.sh $SLAB
  echo "Finished! Type [Enter] to continue."
  read
  exit
fi


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
tau.P=2               ! konstanta barostatu
rescale="ZCM"         ! expand/shring box in the z-direction only

init="start"          ! start from the previously recorded configuration
EOF

cat T.def >> $SLAB.def
cat box.def >> $SLAB.def 
echo "L[0]=x L[1]=y L[2]=3.6*x;" >> $SLAB.def

echo
echo -n "Now the whole box will be thermostatted to temperature "
cat T.def
echo "The x,y-sizes are kept, the box z-size will fluctuate."

echo ">>>>>> simulation parameters:"
cat > $SLAB.get <<EOF
init="start"
dt.plb=1 tau.rho=0 tau.P=2 no=5000;
EOF

cat $SLAB.get

echo
echo "The calculations will be startted in the background."
echo "It is possible to watch the running simulation."
echo "(The stop request will be sent later.)"
echo "[Enter]=continue, [Ctrl-C]=interrupt"
read

rm $SLAB.plb $SLAB.mol
$COOK nacl $SLAB &

echo
echo "Wait a minute or two and start Step A11"

echo "[Enter]=continue"
read
