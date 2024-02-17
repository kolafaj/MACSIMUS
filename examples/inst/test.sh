#!/bin/bash
echo =============================================================================
echo "Test of your MACSIMUS installation."
echo "This is your current environment:"
echo "  BLENDPATH=$BLENDPATH"
echo "  GUI=$GUI"
echo "  PATH=$PATH"
echo "If undefined, try install.sh again and follow the instructions!"
echo
echo "-------------------------------------------------------------"
echo "Now we will run blend to create the force-field file spce.ble."
echo "Type ENTER -- a window with a water molecule will appear."
echo "To continue:"
echo "  - click [finish] in the control panel (rightclick to get help)"
echo "  - or type '.' in the window with the molecule"
echo "If not finished properly, you will get error \"table species not found\" later."
read

function checkinst () {
  type $1 && return
  echo "Bad MACSIMUS installation - '$1' not found!"
  echo "Check whether '$1' has been compiled and is in your path:"
  echo "PATH=$PATH"
  exit
}

[ -e spce.ble ] && rm spce.ble
[ -e SPCE.mol ] && rm SPCE.mol SPCE.gol SPCE.plb

checkinst blend

cat > SPCE.che <<EOF
water SPC/E water model
parameter_set = sea
Hp0.4238 Hp0.4238
    \    /
     On.8476
EOF
blend -g -o spce -h SPCE.che
echo
if ! fgrep "*1-4" spce.ble ; then
  echo "You have not closed 'blend' correctly - incorrect spce.ble was generated!"
  echo "Type ENTER, I'll fix it (by running 'blend -o spce -h SPCE.che')"
  read
  rm SPCE.{mol,gol}
  blend -o spce -h SPCE.che
fi

echo "blend finished"
echo

[ -e spce.def ] && rm spce.def SPCE216.*

# number of particles
N=216

cat > spce.def <<EOF
N[0]=$N         !  number of molecules (quite small for a fast test)
cutoff=9.3       ! [AA] Ewald r-space cutoff
LJcutoff=cutoff  ! [AA] Lennard-jones (or similar) smooth cutoff
rho=1000         ! [kg/m3] density; box will be calculated
equalize.mol=0.8 ! masses of atoms equalized (O=8, H=5), longer step allowed
noint=4          ! MD steps/cycle
h=0.01/noint     ! [ps] MD timestep (4 in 0.1 ps cycle = quite cheap)
dt.plb=0.1       ! [ps] how often to write the trajectory (spce.plb)
dt.prt=0.1       ! [ps] how often to write info to spce.prt
T=298            ! [K] temperature
;
EOF
echo "--- Input file spce.def: ---"
cat spce.def

cat > SPCE$N.get <<EOF
init="crystal" ! start from regular "crystal", random orientations
thermostat="Berendsen" ! thermostat type
tau.T=0.15     ! [ps] thermostat time constant (short to start)
no=200;        ! 100 cycles to fast cool 
thermostat="Nose" ! Nose-Hoover (better thermostat)
tau.T=0.2      ! [ps] thermostat time constant
no=100;        ! 100 cycles to further equilibrate
EOF
echo "--- Input file SPCE$N.get: ---"
cat SPCE$N.get

echo
echo "Type ENTER to start simulation (will finish in < 10 s on a modern computer)"
read

checkinst cookewslc

if ! cookewslc spce SPCE$N ; then
  echo "'cookewslc' failed or has not stopped correctly,"
  echo "try $PWD/$0 again."
  exit
fi

echo
echo "Type ENTER to show the trajectory..."
read

checkinst show

show SPCE$N -I%i
echo
echo "Type ENTER to show selected quantities as the function of time [ps]"
echo "Three graphs will appear:"
echo "  - Etot = kinetic + potential energy [K]"
echo "  - Tkin Tin Ttr = kintic temperature, rotational, translational [K]"
echo "  - P = pressure [Pa]"
echo "Click [killall] or type 'K' in any plot window to kill them all"
read
if ! type showcp ; then
  echo bad MACSIMUS installation - showcp not found!
  echo check whether showcp has been compiled and is in your path
  echo your current PATH=$PATH
  exit
fi

checkinst plot
checkinst showcp

showcp -p SPCE$N Etot Tkin P
echo "Type ENTER to simulate again and calculate radial distributions"
echo "functions and running coordination numbers."
read
cat > SPCE$N.get <<EOF
init="start"      ! start from previous configuration
thermostat="Nose" ! Nose-Hoover (better thermostat)
tau.T=0.2         ! [ps] thermostat time constant
rdf.cutoff=cutoff ! range for RDF
rdf.grid=25       ! [1/AA] grid (25 bins per 1 AA)
no=300;           ! cycles of the productive run
EOF
if ! cookewslc spce SPCE$N ; then
  echo "cookewslc failed or has not stopped correctly,"
  echo "try $PWD/$0 again."
  exit
fi
echo
echo "HINTs:"
echo "- Zoom in by the left mouse button"
echo "- Click by the right mouse button to get curve labels"
rdfg -g -p SPCE$N.rdf
rdfg -c -p SPCE$N.rdf

echo
echo "Erase all files created during running this test (Y/n)?"
read
if [ "$REPLY" == "y" -o "$REPLY" == "Y" ] ; then
  rm SPCE$N.* spce.* SPCE.*
fi
