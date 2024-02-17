/* molecule-molecule interactions */
extern pot1_t intramol;
extern pot2_t *Potential(int i,int j),intermol,atomatom;

#ifdef SLAB
double mol_wall(vector f[],molecule_t *m,vector *rp);
#endif /*# SLAB */
