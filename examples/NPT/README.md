# Simple test of the NPT ensemble

This example is based on paper *Novel barostat implementation for molecular dynamics* by J. **Janek** and J. **Kolafa**.  The Martyna-Tobias-Klein thermostat/barostat is implemented in MACSIMUS using Verlet/leap-frog, SHAKE, and predictors.

## Preparation

Prepare force field from data in `MACSIMUS/blend/data/gases.par` and `N2.che`:

`$ blend -o n2 N2.che`

Assuming that at least 4 threads are available, tell cook to use them.  It does not make sense to use more for such a small system.

`$ export NSLOTS=4`

## Initialization

Run initialization:

`$ cookewslcP1 n2.ble N256NPT1.get`

Check the results:

* Have a look at protocol `N256NPT1.prt`

* Check the timing:<br >
  `$ fgrep speed N256NPT1.prt`                                             

* Check the convergence profiles (three equivalent ways):<br >
  `$ showcp -p N256NPT1.cp Tkin rho P`<br >
  `$ start N256NPT1.cp`<br >
  from Midnight Commander: open file `N256NPT1.cp`

* If not converged, remove the first block in `N256NPT1.get` up to the first `;` and repeat simulation.

## Productive run

When converged, leave only block starting with<br />
`! productive`<br />
on file `N256NPT1.get` and run the simulation again.  It may take ~ 1–2 hours; the simulation for the paper was longer.

`$ cookewslcP1 n2.ble N256NPT1.get`

## Results

The results can be found at the end of the protocol `file N256NPT1.prt`, or using command:

`$ staprt -nV -kNaevnf N256NPT1.sta`<br />
`$ staprt -n'P(tens) [Pa]' -kNaevnf N256NPT1.sta`<br />
`$ staprt -n'enthalpy*' -kNaevnf N256NPT1.sta`

Where the `-k` arguments correspont do the columns printed: N=variable name, a=average, e-error (1σ), v-variance, n=# of data, f=file
