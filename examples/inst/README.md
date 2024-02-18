# Simple test

This test is run directly from the installation script to check the installation.

The following components are tested:
* force field parameters `blend/data/water.par`
* `blend`, force field builder
* `cookewslc`, MD simulation (serial version with Ewald, linked-cell list)
* `showcp`, analyze the convergence profiles
* `plot`, show graphs
* `show`, show the trajectory

Called from the installation script or separately as:
`$ ./test.sh`
