/* support for optimized water models */

extern struct waters_s {
  int iH,iO,iM;                       /* standard indices in *.ble files */
#if defined(WATERPLUS) && WATERPLUS>1
  struct { sitesite_t HH,OH,OO,OM,MM,MH,dummy; } ss; /* LJ +rdf info
                                      (dummy=no LJ+rdf not measured) */
  struct {     double HH,OH,OO,OM,MM,MH; } qq; /* product of charges */
#else
  struct { sitesite_t HH,OH,OO,dummy; } ss; /* LJ +rdf info
                                      (dummy=no LJ+rdf not measured) */
  struct {     double HH,MH,MM; } qq; /* product of charges */
#endif
} waters;
/* NOTE: M=extra charged site for TIP4P
         M=L=lone pair for ST2
         M=O for TIP3P */

pot_t TIP3P,TIP4P,ST2;

#ifdef WATERPLUS
/* WARNING: for one special project only */
#if WATERPLUS==2
pot_t dljd,CO2,acetone,MeOH;
#endif
#if WATERPLUS==3
pot_t TIP3Pion,ionion;
#endif
/* see also "watercut.h" */
#endif
