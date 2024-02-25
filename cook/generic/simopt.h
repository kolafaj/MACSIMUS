#ifndef SIMOPT
#define SIMOPT

/************************************************************************
  cook* compile-time customization
  See also the configuration script  macsimus/cook/configure.sh,  which
  covers most typical cases and sets dependent switches automatically
  ***********************************************************************
  Macros (switches) #defined in this file apply to the whole cook*
  project with the exception of the general software (in directory `gen')

  Enable (uncomment) #defines for your project!
  Some reasonable defaults are already pre-#defined.

  NOTE: statement `#include "simopt.h"' is in "simglob.h".
  "simopt.h" should be then directly #included only to modules that do
  not #include "simglob.h"  (like "erfc.h")
************************************************************************/

#define PROJECT "COOK-your-choice"
/* project name - any text you find useful (just printed in the head) */

#define COOK
#define Lagrange
/* these must be #defined for the cook project */


/************************** BOUNDARY CONDITIONS ***********************/

/* exactly one of FREEBC, NIBC, COULOMB must be #defined */

/*
#define FREEBC
  Free boundary conditions and 1/r electrostatics instead of cubic periodic
  boundary conditions.  Also turns off Ewald (and Ewald-like, like
  cutoff electrostatics).
  Use -ffast-math option when compiled by gcc (because of sqrt())
  Recommended executable name: cookfree
*/

/*
#define NIBC
  Nearest-image boundary conditions: no spherical cutoff, no Ewald.
  Also turns off Ewald (and Ewald-like, like cutoff electrostatics.
  Incompatible with FREEBC or LINKCELL.
  Recommended executable name: cookni
  WARNING: Not tested for a long time.
*/

/*
#define COULOMB -2
  Electrostatics approximation.
  Implies periodic boundary conditions, or periodic in x,y and slab or
  wall in z (formerly "Ewald-like").  Valid values are:
  -3: Gaussian charges on a spring + Ewald, cf. GAUSSIANCHARGES, SPLINES
  -2: The same as COULOMB=-1 and `exact' erfc [see exacterud_sqrt()]
      is used for exceptions instead of splines to avoid imprecise
      results for short separations.  Necessary if there are large
      charges close together (e.g., mimicking a point dipole),
      for surface tension, and simular.
      1-4 interactions are still calculated by splines.
  -1: Ewald summation, r-space part by hyperbolic splines.
   0: MACSIMUS cheap cut-off and smoothed electrostatics.
      Smoother and usually more accurate than COULOMB=3, suitable also
      for large charges close together.
   1: not used
   2: MACSIMUS cut-off electrostatics by quadratic splines (deprecated).
   3: Fennell and Gezelter [JCP 124, 234104 (2006)] cut-and-shift
      potential by quadratic splines.  Not suitable for large charges
      close together
  Recommended executable names: cookew, cookce, cookfg
*/

/*
#define SPLINES -2
  Spline order for electrostatic calculations, with COULOMB<0:
  -2 : hyperbolic splines for erfc/r, good for standard Ewald, bad for
       Gaussian charges
   2 : quadratic splines (not useful, cf. COULOMB=2)
   3 : cubic splines, required for Gaussian charges, also possible for
       standard Ewald but likely less efficient (cache dependent)
*/

/*
#define QQTAB 1
   Separate spline for each charge-charge pair
   For Ewald r-space forces, needed with Gaussian charges:
   (default = function erfc is splined, then multiplied by q*q)
   QQTAB values:
    0 : for models (as BK3) of COS-type, i.e., central atom has a zero charge
        (error if used for models with a nonzero central charge)
    1 : general for models with atom (center) charge + Drude charge
    2 : as above, optimized for both cases, good for a mixture of sites
        with zero and nonzero center charges
*/

/*
#define GOLD
  DO NOT USE
  WARNING: support removed and/or have not been used for a long time
  Simulation on a `golden plate' - conducting surface
  Compatible with FREEBC and COULOMB>=0 (cutoff electrostatics: periodic
    b.c. in x,y)
  Incompatible with LINKCELL
  Implies WALL (SLAB&4)
*/

/*
#define WALL
  DO NOT USE
  Replaced by SLAB&4, NOT TESTED RECENTLY
  1 or 2 walls at z=0 and z=L are added
  GOLD version is automatically also WALL
*/

/*
#define METAL
  metal force fiels w. electron density
  - only direct evaluation of pair forces supported
  - requires one extra pair-pass
  - see sim/metal (hard-wired gold)
*/

/*
#define SHARPCUTOFF
  turn off LJ-like cutoff; needed by METAL
*/

/*
#define POTCONST
  print LJ potential constants, see sim/lj/sitesite.c and sim/setss.c
*/

/*************************** PERFORMANCE **************************/

/*
#define WORM
  Periodic b.c. modifier: Molecules of maximum length 3/2 of the box size
  are allowed.
  If WORM is not specified then the length limit is half the box
  size, i.e., only small molecules can be simulated.
  Use for long molecules like proteins, polymers, etc.
*/

/*
#define PERSUM
  Cutoff longer than half box is allowed so that the pair sum includes
  several periodic images.  Slow!
  Implies the direct sum method (not LINKCELL).
*/

/*
#define CUT
  Additional if's are addded to tests whether site-site distances are
  within the cutoff.  Slightly more efficient if the cutoff is
  substantially less than half the box size (but usually the LINKCELL
  method is even more efficient).
  With LINKCELL only the initializer is affected.
*/

/*
#define LINKCELL 3
  Linked-cell list method.
  Compatible with most COULOMB versions incl. WALL.
  Incompatible with FREEBC or NIBC
  Recommended executable name: cooklc[.exe] (cooklcew, cooklcce)
  Versions:
  0: Forces summed up in vectors which reside physically in the
     particle control structure (linked to the cell-lists).
     After finishing the r-space part they are copied back to contiguous
     arrays. (New since V2.7e)
  1: Forces are (via indirection) summed to the contiguous arrays (the
     particle control structure contains a pointer).
     with PARALLEL==1 this applies to the 1st of pair only
  3: As above also for the 2nd pair (all forces summed directly in the arrays)
     (if not PARALLEL==1 then LINKCELL=1 is the same as 3 because only
     bit 0 is tested)
     Prior cook V2.7d, only this version was available
  (Which version is more efficient depends mainly on the processor and
  cache architecture)
*/

/*
#define PARALLEL 1
  Shared-memory parallel versions for pthreads
  Available versions:
  0: Serial version.
  1: Linked-cell list parallelized like systolic loop, Ewald parallelized.
 -1: As above without barriers (deprecated: starts/stops threads more often)
  2: Two threads, Ewald k-space || r-space
  3: REMOVED
  Variant (not recommended): if SHM is defined in all modules, fork is used instead of threads.
*/

/*
#define SERIAL
  With PARALLEL>0, the parallel code is used, however, all threads are
  executed serially as normal serial subroutines
  For debugging only, see simglob.h
*/

/******************************* EXTENSIONS *******************************/

/*
#define POLAR
  Support of polarizable atoms.
  #undef POLAR = no polarizability
  #define POLAR FLAG:
   0: only scalar linear polarizability supported.
      More features can be obtained by summing the following flags:
   1: repulsive antipolarizability (shell core model by Tosi etal.)
      all terms included
   2: saturated polarizability
   4: use this if there are also uncharged atoms (faster code)
   8: axial polarizability
  16: special (see constrdaa.c and gear2polg.c)
  32: fluctuating charge (FQ4) support
      [old meaning, modifier of POLAR=1, was replaced by "shellrep"]
  64: incl.all intramolecular induced dipole--induced dipole interactions
      (default = 1-2, 1-3 are excluded)
  combination with LINKCELL is in alpha-testing stage
  Recommended names: cookpol etc.
*/

#define SHAKE 1
/*
  enable/disable Shake which appears as method with order=2 (-m2 option) in
  the Lagrangian formalism; Gear+constraint dynamics is still available.

  undefined: Shake not implemented
  1: simple update in sweeps
  2: info on moved sites kept - only bonds of sites that have moved in
     previous step are updated
  negative: as above and a more complicated algorithm taking into account
     the angle between the old and new constraint is used - and some
     tests are added
*/

#define VERLET 3
/*
  Selects formula for kinetic energy (and pressure tensor) with
  Verlet+Shake integration.  Void with Gear+Lagrangian constraints.
  0: SHIFTED, v(t) = [r(t+h)-r(t)]/h
     O(h), harmonic oscillator: Etot+O(h), correct <Ekin>
     (worst total energy => not recommended)
  1: VELOCITY, v(t) = [r(t+h)-r(t-h)]/2h, Ekin = m/2 v(t)^2
     O(h^2), harmonic oscillator: Etot+O(h^2), <Ekin>+O(h^2)
     (equivalent to the velocity Verlet algorithm)
  2: HARMONIC, Ekin = m/2 * [r(t+h)-r(t)]/h * [r(t)-r(t-h)]/h
     O(h^2), harmonic oscillator: constant Etot, <Ekin>+O(h^2)
  3: LEAPFROG, Ekin = m/4 * [(r(t+h)-r(t)]/h)^2 + (r(t)-r(t-h)]/h)^2]
     O(h^2), harmonic oscillator: Etot+O(h^2), correct <Ekin>
     (recommended)
  4: PREDICTED, Ekin = m/2 * v_pred^2, where v_pred is predicted from
     previous [r(t+h-ih)-r(t-ih)]/h
     (only with Nose-Hoover for Verlet+Shake, not recommended)
  5: EQUIPARTITION, Ekin = 5V1^2 + 5V2^2 + 2V1V2,
     where V2 = [r(t+h)-r(t)]/h, V1 = [r(t)-r(t-h)]/h
     (best equipartition for harmonic diatomics)
  (see c/alg/harm.c and c/alg/nose+verlet.c)
*/

/***************************** MEASUREMENTS ****************************/

/*
#define WIDOM
  WARNING: have not been tested recently, WALL removed
  Widom insertion particle method
  there are two versions: for WALL and plain (bulk)
*/

/*
#define ANCHOR
  Allows keeping one site fixed (use get(anchor.i) in data)
  and measures the force on it
*/

/*
#define SLAB 0
  - calculate the z-density profiles, both center-of-mass based and atom based
  - add slab forces
  - add the post-processed slab-based cutoff corection
  - slab bias function available (see center.z1)
#define SLAB 1
  - slab equilibrium in 1 simulation: Lx,Ly as polynomial(T)
#define SLAB 2
  - cleaving
#define SLAB 4
  - LJ wall, replacing the WALL switch, NOT TESTED RECENTLY
#define SLAB 8
  - experimental: testing of "mapped averaging"

NB: efficiency is not affected if these options are included and not used

VERSIONS 2.8b and older:
[[ OLD VERSION, now always present ]]
[[ #define SLAB 1 ]]
[[  - in addition, add slab bias function (see center.z1) ]]
[[ #define SLAB 8 ]]
[[  - changed to SLAB=1 ]]
*/

/*
#define PRESSURETENSOR 3
  include direct pressure (stress) tensor calculations incl. Ewald
  0 = no pressure tensor
  default (PRESSURETENSOR undefined) = 3
  Sum the following flags:
   1: virial (static) part incl. constraint forces, diagonal only
   2: kinetic part, diagonal only
   3: full diagonal (default)
   4: calculate also off-diagonal terms Pxy,Pyz,Pxz"
      (if not set, only diagonal terms Pxx,Pyy,Pzz are calculated)
   8: calculate also center-of-mass based (molecular) tensor
      - applies for the kinetic part only (use with 2)!
      (if not set, only site-based part is calculated)
  16: calculate curtosis: PRESSURETENSOR=27 recommended
  32: debugging of dependant forces only: term is transposed
  64: debugging of dependant forces only: average of standard+transposed
      see norm.c (all versions should be identical)
  32+64: do not use (the same as 32)
  none of 32,64: standard version (recommended)
For direct calculation of surface tension, use 3 for models with
constrained bonds.
1 is good enough for vibrating bonds or monoatomic models in isotropic cube.
*/

/*
#define DIHHIST 4
  Dihedral angle distribution support
   2: only the central bond distinguishes dihedrals
   4: all 4 atoms distinguish dihedrals
  -1: all dihedrals in all molecules different:
      support for selection (SIMNAME.ddh) enabled
   1: one summary histogram for all dihedrals
*/

/*
#define XSECTION
  Enables cross-section calculations.
  See sfdx.c, function Xsection() for details
*/

/*
#define RGYR
  Radius of gyration and end-to-end separation are calculated and
  stored in CP
*/

/*
#define SHEAR
  Enable measuring of shear viscosity (variable shear must be set).
*/

/*
#define CLUSTERS
  cluster analysis, file SIMNAME.cli needed
  (to define connectiveness via distances, etc.)
*/

/*
#define LOG
  In the direct pair versions (not LINKCELL), allows for separate
  recording of Lennard-Jones and similar terms and elecrostatic terms.
  Also, keyword No.first is enabled (some functionality of former option -j).
*/

/************************** SPECIAL VERSIONS *********************/

/*
#define ECC
  Electronic Continuum Correction
  experimental code
*/

/*
#define HARMONICS 3
  Harmonic expansion coefficients + full correlation function
  3=water, 2=linear molecule, (1=methanol)
  Special for water etc. models, see:
     sim/harmonic.h
     sim/harmwat.c, harmco2.c ...
     util/harmg.c, harmgco2...
*/

/*
#define STARS
  gravitational forces instead of Coulomb (sign changed), FREEBC required
*/

/*
#define GAUSSIAN
  calls Gaussian to calculate force and energy, FREEBC required
  many measurements are not available
  atom names must be understood by Gaussian input - might require ble-file edit
  see the manual
*/

/*
#define MARK
  special version for energy partitioning (e.g., residue-residue)
  not tested recently
*/

/*
#define BJERRUM
  not public
  module minimize must be added to metamake
*/

/****************************** SYSTEM ****************************/

/*
#define FLOAT
  The configuration, potential tables, and the potential and forces
  calculations are in single precision.
  Deprecated!
  Generally not recommended unless very simple systems.
*/

/*
#define IFLOAT
  Use `#define IFLOAT' if you want to have only intramolecular tables in
  single precision.
  Deprecated!
  FLOAT implies automatically IFLOAT
*/

/*
#define REDUNITS
  Reduced units - all unit factors are 1.
  Never used.
*/

/*
#define CLOCK
  Additional code to measure time using function clock()
  in the way compatible with the FMM project
  System-dependent code!
*/

/* for obsolete and special switches, see simopt-o.h */

#endif /* SIMOPT not defined */
