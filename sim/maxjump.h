/* 
   to interface diffusivity and conductivity measurements in modules
   sim/simmeasx.c (former sim/sfdx.c) and util/plb2diff.c 
*/
struct maxjump_s {
  double xi[3]; /* max jump, separately x,y,z */
  int frame[3]; /* frame */
  int n[3],i[3];/* molecule,site */
  int no;       /* f,init removed in V3.0a */ 
};
