/* 
   To interface diffusivity and conductivity measurements in modules
   sim/simmeasx.c (former sim/sfdx.c) and util/plb2diff.c.
   NB: "jump" here means the displacement extended by following
       the trajectory over periodic boundaries
*/
struct maxjump_s {
  vector r;     /* (max) jump (displacement), separately x,y,z, with sign */
  int frame[3]; /* frame of plb-file */
  int n[3];     /* molecule number */
  int to;       /* max. frame (reread.to) */
  /* f,init removed in V3.0a */ 
};
