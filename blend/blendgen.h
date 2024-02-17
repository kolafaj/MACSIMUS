/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BLENDGEN.H

This module assigns parameters of atoms, bonds, angles, dihedrals and
impropers to the molecules and generates input file for `cook' and `cooks'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

extern int nsites;   /* total # of different sites used */
extern int clust;    /* cluster # while splitting molecule */

extern site_t *site; /* global image of = spec->site */

#define T(I) site[I].type

/* NOTE: loopnbr moved to blendmed.h */

extern int prtatomoffset;
void prtatom(int i); /* prt atom name */

void listofsites(void); /* list of all site types used */

void build(species_t *spec); /* builds molecule
				(force field, minimize, print .ble) */
