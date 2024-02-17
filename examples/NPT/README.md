# Simple test of the NPT ensemble

This example based on paper *Novel barostat implementation for molecular dynamics* by J. Janek and J. Kolafa: The Martyna-Tobias-Klein thermostat/barostat implemented in MACSIMUS using Verlet/leap-frog, SHAKE, and predictors

Prepare force field from data in `MACSIMUS/blend/data/gases.par` and `N2.che`:

`$ blend -o n2 N2.che`

Assuming that at least 4 threads are available, tell cook to use them.  It does not mae sense to use more for such a small system.

`$ export NSLOTS=4`

Run initialization:

`$ cookewslcP1 n2.ble N256NPT1.get`

Check the convergence profiles:

* Command:

  `$ showcp -p N256NPT1.cp Tkin rho P`

* The same if `start` is installed:

  start N256NPT1.cp

* Or from Midnight Commander:

  open file `N256NPT1.cp`

Check the timing by:

`$ fgrep speed N256NPT1.prt`                                             

*If not converged, remove the first block in `N256NPT1.get` up to the first `;` and repeat simulation.*

If converged, run (~ 1-2 hours; the simulation for the paper was longer):
  `$ cookewslcP1 n2.ble N256NPT1.get`

The results can be found at the end of the protocol file N256NPT1.prt, or using command:

`$ staprt -nV -knaevnf N256NPT1.sta`
`$ staprt -n'P(tens) [Pa]' -knaevnf N256NPT1.sta`
`$ staprt -n'enthalpy*' -knaevnf N256NPT1.sta`

Where the `-k` arguments mean: n=name, a=average, e-error (1σ), v-variance, f=file
