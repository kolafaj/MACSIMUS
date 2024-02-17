#!/bin/bash

source setup.sh ce

[ x$SIZE == x ] && SIZE=5
X=`ev "$SIZE*128+185"`
Y=`ev "$SIZE*120"`

echo "%%%%%%%%%%% STEP B02 - equilibrate in the NPT ensemble %%%%%%%%%%%"

if [ ! -e solution.cfg ] ; then
  echo
  echo "ERROR! File solution.cfg is missing, repeat Step B01.sh"
  sleep 3
  exit
fi

cat > solution.get <<EOF
init="start"
corr=3 ! cutoff corrections + kinetic pressure
tau.P=2 ! Berendsen barostat
no=500*a; ! = 5 ps for simple solute, 10 ps for endofullerenes
EOF

echo ">>>>>> simulation parameters:"
cat solution.get
rm solution.cp
echo
echo "Let's simulate - type [Enter]"
read

[ "$COOK" == "" ] && COOK=cookce

while [ x$R == x ] ; do

  $COOK ff solution

  if [ ! -e solution.cp ] ; then
    echo "ERROR! File solution.cp is missing."
    echo "Type [Enter]..."
    read
    exit
  fi

  echo
  echo "Now I will show the time dependence of temperature, energy, and density."
  echo "If the graphs do not converge, repeat equilibration!"
  echo "[Enter]=continue, [Ctrl-C]=interrupt"
  read
  nice showcp -p5 solution Tkin Epot rho

  echo
  if ./yes.sh "Repeat equilibrating with start=\"append\" (Y/n)?" y ; then R="" ; else R=n ; fi

cat > solution.get <<EOF
init="append"
tau.P=2   ! barostat time constant
no=500*a  ! the number of steps by 0.01 ps
;  
EOF

done
