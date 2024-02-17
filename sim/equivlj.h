/* 
  To be #included by sim/XXX/sitesite.h to struct pairaux_s if SLAB.
  The parameters refer to the equivalent LJ potential (the same minimum)
  as a non-LJ potential in sim/XXX.
  The LJ potential is assumed in the form
    uLJ(r) = 4*LJeps*[(LJsig/r)^12 - (LJsig/r)^6]
  and Sq=LJsig^2, E4=4*LJeps.
*/
real LJeps,LJsig,Sq,E4;
