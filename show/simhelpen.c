const char *HELP="\
This program shows transitions between solid, liquid, and gas states\n\
as well as other molecular phenomena, as modeled by computer simulation.\n\
\n\
Molecules attract each other when they are close together.\n\
In our model this happens when the atom-circles almost touch.\n\
If the molecules are even closer, they begin to repel.\n\
In our model the molecules repel when the atom-circles overlap.\n\
If the atoms are far away, their attraction becomes negligible.\n\
\n\
Try to cool or warm up the system - it may take a while!\n\
The solid phase (crystal) is stable around temperature of 0.05,\n\
liquid around 0.15, gas around 1.\n\
\n\
The temperature (in \"arbitrary units\", not degrees):\n\
  T = thermostat temperature\n\
  Tk/Td = actual temperature (fluctuates during the simulation)\n\
\n\
The graph shows the radial distribution function which represents\n\
the probability of finding a pair of molecules a goven distance apart,\n\
normalized so that for interactionless particles (ideal gas) it is 1.\n\
\n\
----------------------------------------------------------\n\
* Summary of used acronyms and more hot keys: press h\n\
* In case of flickering: type a\n\
* Quit simolant: ESC (type twice)\n\
* Get help for given function button by RIGHT mouse click.\n\
----------------------------------------------------------\n\
\n\
Now continue by any key...";

const char *SPEEDHELP="\
These buttons control the show and simulation speed.\n\
Hot keys: s S\n\
\n\
In detail:\n\
- r=+0 means that after one MD step or one attempted move per\n\
  particle the results (radial distribution function, pressure) are\n\
  calculated and the configuration shown\n\
- r<0 means that there are delays included\n\
- r>0 means that not all steps are shown\n\
  and the results are not calculated so often";

const char *THELP="\
Press [cool] or move the wheel to decrease temperature.\n\
Press [warm] or move the wheel to increase temperature.\n\
Hot keys: t T\n\
\n\
The temperature set is the thermostat temperature.\n\
It takes a while for the system until it reaches the\n\
selected temperature.\n\
\n\
For experts:\n\
Temperature is a parameter of the Monte Carlo method.\n\
Imagine that each molecule touches a thermostat while moving.\n\
High-energy positions are more probable at high temperatures,\n\
low-energy positions are stable at low temperatures.\n\
\n\
In the molecular dynamics, the total energy (kinetic+potential)\n\
is conserved.  To keep the temperature constant, all velocities\n\
are changed a bit at every step to approach the selected temperature.\n\
\n\
Note: hot key e toggles the constant energy/constant temperature\n\
simulation\n";

const char *WALLHELP="\
By pressing these buttons you can switch between\n\
attractive and repulsive walls.\n\
Hot keys: x X y Y\n\
\n\
Attractive walls are green, repulsive walls are red.\n\
The wall interaction parameter can be set during start.";

const char *NBRHELP="\
A coordination number is shown for each molecule.\n\
Hot key: n\n\
\n\
The coordination number is defined as the number of neighbors\n\
upto the center-center distance of 1.5 of molecule diameter.\n\
\n\
In the hexagonal 2D crystal the coordination number is 6\n\
(more in a real 3D crystal, according to the crystal lattice).\n\
In the 2D liquid it is around 5 (10-12 in a 3D liquid),\n\
and in gas it is small.\n\
\n\
The simulation restarts in a few seconds...\n\
\n\
For experts:\n\
This number can be calculated from the radial distribution\n\
function by integration to distance 1.5 over element rho*2*pi*r*dr.";

const char *METHODHELP="\
Two simulation methods are available:\n\
\n\
Monte Carlo samples all possible positions of molecules using\n\
random moves.  A random displacement is accepted or rejected\n\
according to the energy of the new position and temperature.\n\
E.g., it is rejected should the molecules overlap.\n\
\n\
Molecular dynamics calculates trajectories of molecules in time by\n\
numerical integration of the Newtonian equations of motion.\n\
Molecules move \"as in the real world\" - just many times slower.\n\
\n\
The methods are toggled by button [MC/MD] (or hot key m).\n\
At start [auto] (hot key u) is selected which means that\n\
a suitable method is chosen automatically.\n\
\n\
Note: the displacement length in MC and the integration\n\
step in MD are adjusted automatically.  This can be\n\
changed by hot keys d D.";

const char *GASHELP="\
By this button you will make gas at temperature T=1.\n\
Hot key: F5\n\
Try to refrigerate and observe the changes!";

const char *LIQUIDHELP="\
This button produces a small droplet.\n\
It will solidify on freezing and evaporate on heating up.\n\
Hot key: F6";

const char *CAPILHELP="\
Capillary action under a microscope.\n\
Hot key: F7\n\
\n\
The side and bottom walls are made attractive and temperature\n\
is set to T=0.15 which is in the liquid range.\n\
In addition, small gravity is added.\n\
\n\
More: \n\
- gravity can be changed by hot keys g G\n\
- attractive/repulsive walls are toggled by x X y Y\n\
\n\
Hint: falling drop:\n\
- set the temperature to about 0.13\n\
- make the upper wall attractive, other walls repulsive\n\
- use positive gravity until all molecules stick to the upper wall\n\
- set gravity to -0.01: a droplet will fall";

const char *IDHELP="\
Ideal crystal.\n\
Hot key: F9\n\
Our 2D model crystallizes in a hexagonal lattice.\n\
Try to melt it - type [warm] several times.";

const char *EDGEHELP="\
Ideal crystal with a (2D) dislocation: a part of one row of\n\
atoms is missing.  The dislocation travels during the simulation\n\
and eventually disappears at the border of the crystal.\n\
Hot key: F10\n\
\n\
If you do not see clearly the defect, use \"show: [neighbors]\".";

const char *VACANCYHELP="\
Ideal crystal with a point defect: a vacancy.\n\
The defect travels during the simulation and finally\n\
disappears at the border of the crystal.\n\
Hot key: F11\n\
\n\
If you do not see clearly the defect, use \"show: [neighbors]\".";

const char *INTERHELP="\
Ideal crystal with a point defect: an interstitial atom.\n\
The defect travels during the simulation and finally disappears\n\
at the border of the crystal.\n\
Hot key: F12\n\
If you do not see clearly the defect, use \"show: [neighbors]\".";

const char *GHELP="\
Press [] to increase gravity (also hot key G).\n\
Press [] to decrease gravity (also hot key g).\n\
Press [0] to set gravity to zero (also hot key z).\n\
\n\
For gas (T>1) and gravity set you will see that the gas \n\
is denser at ground according to the barometric formula.";

