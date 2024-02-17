extern int ble_file;
int findsite(species_t *spec,char *id,int chk);
void moledit(species_t *spec);
void react(species_t *spec,char *REAname);
void autobonds(species_t *spec,int rg);

#ifndef TINY
void readjet(species_t *spec,bond_t **b0ptr,fixdih_t **df0ptr);
#endif
