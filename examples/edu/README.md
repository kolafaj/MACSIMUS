# Three educational projects

All projects are designed as a set of scripts which guide a student through all needed steps in a simulation experiment. Knowledge of input files syntax or even linux (if controlled by the Midnight Commander) is no required.

See `simenw2.pdf` included for more details.

## Setup
If you have several processors (not cores), set the environment variable:
`$ export NSLOTS=2`
The scripts will use a parallel version of cook.
You may check the number of processors you have by:
`$ less /proc/cpuinfo`

Enviroment variable SIZE determines the sizes of windows.
The default is SIZE=4. If the graphs are too small or large, change SIZE (example):
`$ export SIZE=6`

Midnight commander (MC) with MACSIMUS associated extensions is recommended.

## *Task A* Melting point of NaCl

Melting point of a model of NaCl by simulating direct coexistence crystal/melt in the slab geometry.

Run scripts:
`$ A01-prepare-Na4Cl4.sh`
`$ A02-replicate.sh`
...

Note: This task is designed for a group of students working in parallel in subdirectories, each student running the simulations at different temperature.  Then, `A12-show-all-profiles.sh` shows all obtained convergence profiles (density vs. time).

## *Task B* Structure of water around a simple solute

A student chooses a solute, runs simulations, analyses radial distribution functions and running coordination numbers, and watches the first/second hydration shell incl. hyddrogen bonds highlighted.

Run scripts:
`$ B01-NVT-start.sh`
...

## *Task C* Coalescence of droplets.

Two small droplets are created and let coelesce. A student analyses the change of temperature caused by the surface energy.

Run scripts:
`$ C01-box-of-water.sh`
...
