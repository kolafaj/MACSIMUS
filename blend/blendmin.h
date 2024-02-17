extern bond_t *b0;
extern angle_t *a0;
extern torsion_t *d0,*i0,*ar0,*cis0;
double minimize(species_t *spec,int print,enum keep_e mask);
#ifndef TINY
extern fixdih_t *fd0;
extern int constraints;
extern int aty_n;
void anglemsd(species_t *spec);
void normalmodes(species_t *spec);
void essential(species_t *spec);
void inertiamatrix(char *IG,species_t *spec);
void virial(species_t *spec1,species_t *spec2,species_t *spec);
#endif
