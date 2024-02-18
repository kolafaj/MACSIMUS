# Brine in the slab geometry using the polarizable BK3/MAH model

Motivated by the work with P. Jedlovszky, see [this paper](https://doi.org/10.1021/acs.jpcb.3c02641) and references therein.

The BK3 water model and compatible salts are based on TIP4P water geometry and Gaussian Drude charges.  The models are polarizable. See [this paper](http://doi.org/10.1063/1.4895129) and references.

## Compiling `cook*`

To configure the `cook*` version, run script `configure.sh` from `MACSIMUS/cook.`; e.g.:

`$ ./configure.sh gcpol/P1`

During configuration, choose defaults except for:

===== Select algorithm for ELECTROSTATIC calculations:<br >
-3

===== Select POLARIZABILITY support:<br >
0

===== Select DETAILS of boundary conditions:<br >
s

===== The NON-BONDED forces are:<br >
busing

===== Are all your molecules smaller than half the box size (y/N)?<br >
y

===== Select PARALLELIZATION (using linux threads):<br >
1

You should obtain `cookgcp0slcP1` (gc=Gaussian Charges, p0=POLAR=0, s=SLAB, lc=LINKLIST, P1=PARALLEL=1).

Also, check the number of available threads on your computer but do not use more than 4 for this system:

`$ export NSLOTS=4`

## Force field

Bug: There is no `.par` file with the  Kiss-Baranyai force field.
Nevertheless `NaCl.ble` is available.
Watch the order of species (keyword `species`).

Also, mol-files for both ions and water are prepared.

## Slab

Run the initialization:

`$ cp 1-slab.get slab.get`<br >
`$ cookgcp0slcP1 NaCl.ble slab`

This takes 7 minutes on Lenovo T480 with i5-8250U CPU @ 1.60GHz.
When done, check the configuration (three options):

* `$ show -Ii slab.plb` # -Ii automatically starts playing
* `$ start slab.plb` # if `start` is installed
* click `slab.plb` from Midnight Commander

Also, check the convergence (three options):

* `$ showcp -p slab Tkin rho P pmax`
* `$ start slab.cp` # if `start` is installed
* click `slab.cp` from Midnight Commander

## Optimizing ASPC

MACSIMUS uses [ASPC](https://doi.org/10.1002/jcc.10385). The default setup guarantees stability, but better precision can be obtained by adjusting the parameters; at the same time, equilibration will proceed.

* `scf.omega` = mixing (damping) parameter. The higher, the more accurate induced dipole, but lower stability. It should be slightly smaller than the stability limit.

* `scf.eps` = accuracy before a second iteration is calculated.  Should be as small as possible so that only one iteration is needed in most cases and only rarely a second one is calculated.

First, use convergence profile of variable `pmax` by 1 step (not blocked):

`$ showcp -p -b1 slab pmax`

It is the maximum error of one iteration of induced dipoles; its monitoring was requested by keyword `pmax` in file `slab.cpi`.  It will fluctuate wildly at start, but then it should converge.  Determine its maximum from about the the second half of the convergence profile.  Write 2x or 3x this value as `scf.eps` to file `slab.def`. Then run:

`$ cp 2-slab.get slab.get`<br />
`$ cookgcp0slcP1 NaCl.ble slab`

Then, write the number reported to `slab.def`.

## Equilibrating

`$ cp 2-slab.get slab.get`<br >
`$ cookgcp0slcP1 NaCl.ble slab`

Watch `pmax` and variables from `slab.prt` named `polar*`. Adjust `scf.eps` and `scf.omega` if needed. If there are peaks on `pmax` over `scf.eps`, decrease `scf.omega`.  Repeat Equilibrating if needed.

The convergence profile of `Etot` includes the Nose extended energies and may have a trend because the method is not exactly time-reversible.

## Productive run

Set variable `no` in `slab.get` and run.

## Analysis

* The calculated z-profiles of number density and charge have extension `.z`. These are ASCII files; to show them using MACSIMUS plot:

  `plot slab.cm.AA-3.z`

  Use button [+] next to `col:` or hot key `]` to change the column shown, [recalc] or `F9` to refresh scaling.
  
* For the surface excess formula, see [Equation (2)](https://doi.org/10.1021/acs.jpcb.0c05547).
* In `slab.prt` you will find the surface tension.
