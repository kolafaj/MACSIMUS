void printdefines(FILE *f)
{
#ifdef NONBOND
  fprintf(f,NONBOND);
#endif
#ifdef REP12
    fprintf(f,REP12);
#endif
#ifdef METAL
    fprintf(f," METAL");
#endif
#ifdef VERLET
    fprintf(f," VERLET=%d",VERLET);
#endif
#ifdef SHAKE
    fprintf(f," SHAKE=%d",SHAKE);
#endif
#ifdef COULOMB
    fprintf(f," COULOMB=%d",COULOMB);
#endif
#ifdef NIBC
    fprintf(f," NIBC");
#endif
#ifdef FREEBC
    fprintf(f," FREEBC");
#endif
#ifdef PRESSURETENSOR
    fprintf(f," PRESSURETENSOR=%d",PRESSURETENSOR);
#endif
#ifdef POLAR
    fprintf(f," POLAR=%d",POLAR);
#endif
#ifdef STARS
    fprintf(f," STARS");
#endif
#ifdef GAUSSIANCHARGES
    fprintf(f," GAUSSIANCHARGES");
#endif
#ifdef QQTAB
    fprintf(f," QQTAB=%d",QQTAB);
#endif
#ifdef SPLINE
    fprintf(f," SPLINE=%d",SPLINE);
#endif
#ifdef LINKCELL
    fprintf(f," LINKCELL=%d",LINKCELL);
#endif
#ifdef PARALLEL
    fprintf(f," PARALLEL=%d",PARALLEL);
#endif
#ifdef SLAB
    fprintf(f," SLAB=%d",SLAB);
#endif
#ifdef SLIT
    fprintf(f," SLIT");
#endif
#ifdef CUT
    fprintf(f," CUT");
#endif
#ifdef WORM
    fprintf(f," WORM");
#endif
#ifdef PERSUM
    fprintf(f," PERSUM");
#endif
#ifdef GOLD
    fprintf(f," GOLD");
#endif
#ifdef WIDOM
    fprintf(f," WIDOM");
#endif
#ifdef ANCHOR
    fprintf(f," ANCHOR");
#endif
#ifdef RGYR
    fprintf(f," RGYR");
#endif
#ifdef XSECTION
    fprintf(f," XSECTION");
#endif
#ifdef SHEAR
    fprintf(f," SHEAR");
#endif
#ifdef CLUSTERS
    fprintf(f," CLUSTERS");
#endif
#ifdef HARMONICS
    fprintf(f," HARMONICS=%d",HARMONICS);
#endif
#ifdef DIHHIST
    fprintf(f," DIHHIST=%d",DIHHIST);
#endif
#ifdef PERSUM
    fprintf(f," PERSUM");
#endif
#ifdef LOG
    fprintf(f," LOG");
#endif
#ifdef BJERRUM
    fprintf(f," BJERRUM");
#endif
#ifdef ECC
    fprintf(f," ECC");
#endif
#ifdef MORSE
    fprintf(f," MORSE (warning: expanded)");
#endif /*# MORSE */
}
