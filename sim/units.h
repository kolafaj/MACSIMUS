#ifndef UNITS_H
#  define UNITS_H

/*
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                       units of measurement                      %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MACSIMUS may records the results either in the internal program units or
SI units

In many cases the user output is in SI units, which is indicated in
[brackets].  If no units are given, the internal program units are used.

WARNING: The author of this program thinks that anybody using non-SI units
should be fined and if does not pay then hanged by legs in a chilly draft.
Unfortunately, direct using of SI units would cause over- and underflows in
calculations so that special units had to be used internally in the program.
In addition, some multiplications are spared if the unit of charge is chosen
so that the electrostatic energy is charge*charge/length

The values below Values based on 05/2019 SI definitions and the CODATA 2018
value of the fine-structure constant:
  alpha=1/137.035999084(21).  
Hence, eps0=8.8541878128(13)e-12 F/m.
(Note that https://doi.org/10.1038/s41586-020-2964-7 (2020) gives
  alpha=1/137.035999206(11)
differing by 5 sigma)

===============================================================================
quantity    | program unit |value in SI and other units
------------+--------------+---------------------------------------------------
temperature |K             |1 K
length      |AA            |AA = Angstrom = 1e-10 m
time        |ps            |1e-12 s
energy      |k*K           |1.380649e-23 [J] 
                           |  = 8.31446261815324 [J/mol]
            |              |  = 1.987204258640832 [cal/mol K]
velocity    |AA/ps         |100 m/s
mass        |k*K/(AA/ps)^2 |1.380649e-27 [kg] = 8.31446261815324 [g/mol]
mass density|k*K*ps^2/AA^5 |1380.649 [kg/m3]
pressure    |k*K/AA^3      |13806490 [Pa]
charge      |\(4*pi*eps0*k*K*AA)     
            |              |SI: 3.919412183482151e-22 [C]
            |              |CGS: 1.175010212721575e-12 [cm1.5 g0.5 s-1]
dipole      |AA*\(4*pi*eps0*k*K*AA)
moment      |              |SI: 3.919412183482151e-32 [C m]  
            |              |CGS: 0.01175010212721575 [D]
===============================================================================

NOTES:
  \ = sqrt
  k = Boltzmann constant = 1.380649e-23 J/K
  (loosely, energy is measured in Kelvins)
  AA = Angstrom
  1 cal = thermochemical calorie = 4.184 J
  program electrostatic energy is q1*q2/r
  charge of electron is 408.779826922048 charge units, e^2=167100.9468984195 charge units^2
  1 mol = 6.02214076e23 (Avogadro number)
  SI polarizability volume = SI polarizability/(4 PI epsilon_0)
  The values are consistent within error bars
  1/(4 PI eps_0) = {c}^2 1e-7 = 8987551792.3(13) [Jm/C^2=m/F]
  1 e = 1.602176634e-19 C = 96485.33212331 C/mol (Faraday constant)
  e^2/(4 PI eps_0)/AA = 2.30707755234(34)e-18 J = 1389354.57644(21) J/mol 
    = 332.06371330(5) kcal/mol
  1 D = sqrt(2e-49[N m4]*pi*e^2/(alpha*h*c)) = 3.33564095107(26)e-30 [C m]
  1 e AA = 4.80320471388(37) D
  1 C m = 2.99792458082(23) D

(See also evu and .evudata)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

/*** basic SI factors ***/
#    define massunit 1.380649e-27 /* kg */
#    define lengthunit 1e-10       /* m */
#    define timeunit   1e-12       /* s */
#    define chargeunit 3.919412183482151e-22 /* C (watch 4 pi eps0) */
#    define energyunit 1.380649e-23 /* J */
#    define Avogadro (1e-3*munit/massunit) /* 1/mol */
/* Examples:
   1) viscosity has dimensionality  mass length^-1 time^-1 = pressure*time
      the program unit of viscosity is thus
        viscunit=massunit/lengthunit/timeunit
                =Punit*timeunit
                =1.380649e-05 Pa.s
      Example: from cook*, we get viscosity=123 [p.u.]
        viscosity = 123*1.380649e-05 Pa.s
        The same calculation using evu: 123[-Pa.s]
        
   2) electric current density has dimensionality  charge length^-2 time^-1
      the program unit of current density is thus
        Junit=chargeunit/lengthunit^2/timeunit
             =3.919412183482151e+10 A/m2
      Example: from cook*, we get Jx=0.1 [p.u.]
        Jx=0.1*3.919412183482151e+10 A/m2
        The same calculation using evu: 0.1[-A/m2]*/

/*** some derived program units in SI ***/
#    define Munit     0.831446261815324 /* g/mol */
#    define Eunit     8.31446261815324  /* J/mol */
#    define Punit     13806490.         /* Pa */
#    define rhounit   1380.649          /* kg/m^3 */
#    define forceunit 1.380649e-13      /* N */
#    define Debye 0.01175010212721575 /* this is 1 Debye in program units */
#    define kcal  4184.0              /* J */
#    define electron 408.779826922048 /* charge of electron in program units */

/* h/k for normal mode frequencies: 
  h/k=4.79923734e-11 [s K] (Codata 2010)
  h/k=4.799244662e-11 [s K] (Codata 2014) 
  05/2019: [s K] */
#    define PLANCK_BOLTZMANN 4.799243073366221e-11

/* example:: 
   char info[]="p.u. of dipole moment = " STRING(Debye) " D"; */
#    define STRING(X) #X

#endif /*# UNITS_H */
