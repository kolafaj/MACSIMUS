/*
  %%%%%%%%%%%%%%%%%%%%%%% reading/writing rsd/mol/cfg-files %%%%%%%%%%%%%%%%%%%
*/

extern enum moltype_e test_cterp,test_nterp,test_patch;
/* reads name.rsd file */
site_t *readRSD(char *name,enum moltype_e *moltype_p,int resno,char chain); 

void writeMOL(void); /* Writes the molecule file *.mol */

void write3D(char *mode,int opt); /* writes the configuration (text or bin) */
extern int retonerror; /* returns 0 on some read errors */ 
int read3D(char *mode,int opt,int no); /* reads the configuration */
