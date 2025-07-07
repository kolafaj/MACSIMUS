/*
  measuring: CP, RDF and dihedral angle distribution (#ifdef DIHHIST)
  dipole moments, radius of gyration A
*/

#include <time.h>
#include "ground.h"
#include "sds.h"
#include "varfile.h"
#include "statics.h"
#include "simglob.h"
#include "simmeas.h"
#include "simdef.h"
#include "simils.h"
#include "pakcp.h"
#include "units.h"
#include "forces.h"
#include "interpot.h"
#include "setss.h" // WIDOM only

#ifndef SS_MEASURE_rep
#  define SS_MEASURE_rep /**/
#endif /*# SS_MEASURE_rep */

/*
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                              convergence profile
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#include "cpmark.h"

int NCP;      /* # of items in a record - to grow */

/* OLD VERSION: typedef struct { float Etot,T,Epot,Eintra,P; } CPtype; */

float *CPbuf; /* convergence profile is buffered and stored at once when
                 also .sta and .cfg are stored.  This is useful for
                 recovering the simulations in the case of crash */

/*
  List of variables always recorded in SIMNAME.cp
  This info MUST match with the series of recordCP in maincps.c: watch all #if's!
  Each item must be 4 chars long (or it must be padded by ' ' to 4 chars)
  1st 2 items (Etot and T) cannot be changed and are NOT included here
  Those marked by ? will be modified later.
  Cf. file SIMNAME.cpi
*/
static char CPinfo0[]=
/* col3col4col5col6col7 ; NB: ?M/J now also for FREE */
  "Epot?   ?   Tin Ttr ?M/J"
#ifdef SHEAR
  "Cvel"
#endif /*# SHEAR */
#ifdef POLAR
  "self"
#endif /*# POLAR */
#ifdef DIHHIST
  "gauc"
#endif /*# DIHHIST */
#ifdef RGYR
  "Rgyr"
  "ende"
#endif /*# RGYR */
#ifdef XSECTION
  "Xsec"
#endif /*# XSECTION */
#ifdef WIDOM
  "Widm"
#endif /*# WIDOM */
  ;
static char *CPinfo=NULL;

static int CPnbit;

struct ssd_s *ssd0;

static void prtCPinfo(char *info) /******************************* prtCPinfo */
{
  int i;

  prt_("  1=A=\"Etot\"   2=B=\"T   \"");
  loop (i,0,NCP-2) {
    prt_("%c%3d=%c=\"%c%c%c%c\"",
         " \n"[i%5==3],i+3,i+3<27?i+'C':'#',
         info[4*i],info[4*i+1],info[4*i+2],info[4*i+3]); }
  _n
}

/*
  names to be listed in SIMNAME.cpi
*/
static struct CPItab_s {
  char *name;  /* name referring to .cpi file */
  char *info;  /* one line of human-readable info */
  double *var; /* pointer to a double variable with the value */
} CPItab[]={
  /* {"KEY","INFO",&double}, */
  /* NEW in V3.3s,V3.4a,V3.6g:
     If "INFO" contains "[J/mol]", the value is converted from p.u.=K to J/mol
     If "INFO" contains "[V/m]", the value is converted from p.u. to V/m
     If "INFO" contains "[Pa]", the value is converted from p.u. to Pa
     (no such automatic mechanism applies to other units)  */
  {"Epnc","potential energy without cutoff corrections [J/mol]",&En.pot},
  {"Eusr","user potential energy (see userforces.c) [J/mol]",&En.usr}, // NEW in V2.9f
  {"elst","electrostatic energy (POLAR: excl. self-energy) [J/mol]",&En.el},
  {"bond","bonded energy [J/mol]",&En.bonded},
  {"Eext","extended degrees of freedom energy (Nose-Hoover, MTK barostat) [J/mol]",&En.ext},
  {"Ekin","kinetic energy [J/mol]",&En.kin},
  {"U",   "internal energy = Epot + Ekin with cutoff corrections [J/mol]",&En.U},
  {"Unc", "internal energy = Epot + Ekin without cutoff corrections [J/mol]",&En.Unc},
  {"H",   "enthalpy = Epot + Ekin + pV with cutoff corrections [J/mol]",&En.H},
  //REMOVED  {"Hnc", "Enthalpy = Epot + Ekin + pV, without cutoff corrections [J/mol]",&En.Hnc},
  {"fix", "potential energy of forces fixing sites in place (-k) [J/mol]",&En.fix},

#ifdef ECC
#error TO BE UPDATED  
  {"Pvir","virial pressure except epsf(V) dependence [Pa]",&En.ECC_Pvir},
  {"PECC","ECC pressure, the same as ECC P [Pa]",&En.P},
  {"Psc","conventional pressure w/o any epsf-based ECC terms [Pa]",&En.ECC_Pscaled},
  {"Pcor","ECC correction (to be added to Pvir) [Pa]",&En.ECC_Pcorr},
#else /*# ECC   */
  {"Pref","pressure as defined by variables virial, corrected [Pa]",&En.Pref},
#endif   /*#!ECC   */
  {"PdV","pressure by virtual volume change (also PdVm/PdVa) w. cutoff corr. [Pa]",&En.PdV.c},
  {"Pevc","pressure (based on el. virial) w. LJ cutoff corr. [Pa]",&En.Pevir.c},
  {"Pevn","pressure (based on el. virial) w/o LJ cutoff corr. [Pa]",&En.Pevir.n},
  {"V"   ,"volume [p.u.=AA^3]",&box.V},

#ifdef LOG
  {"Ep0", "potential energy first No.first molecules [J/mol]",&En.pot0},
  {"EpX", "potential energy  No.first vs. rest molecules [J/mol]",&En.potX},
  {"Ein", "intramolecular potential energy [J/mol]",&En.intra}, // NEW in V2.8c
  {"el0", "electrostatic energy first No.first molecules (w/o self) [J/mol]",&En.el0},
  {"elX", "electrostatic energy No.first vs. rest molecules (w/o self) [J/mol]",&En.elX},
  {"LJ",  "non-bonded (LJ) energy excl. elst. [J/mol]",&En.LJ},
  {"LJ0", "non-bonded energy excl. elst. first No.first molecules [J/mol]",&En.LJ0},
  {"LJX", "non-bonded energy  No.first vs. rest molecules [J/mol]",&En.LJX},
  {"bon0","bonded energy first No.first molecules [J/mol]",&En.bonded0},
#endif /*# LOG */
  /* the reference point = (0,0,0) for FREEBC, otherwise box.center (sf. SLAB) */
  {"CMdr","corrected center-of-mass drift [AA=p.u.]",&En.CMdrift},
  {"vdr", "corrected velocity drift [p.u.]",&En.vdrift},
  {"AMdr","corrected angular momentum drift [p.u.]",&En.AMdrift},
  {"omdr","corrected angular velocity drift [p.u.]",&En.omegadrift},
  {"CMx", "center-of-mass [AA]",En.CM},
  {"CMy", "center-of-mass [AA]",En.CM+1},
  {"CMz", "center-of-mass [AA]",En.CM+2},
  {"LMx", "linear momentum [kg m/s]",En.LM},
  {"LMy", "linear momentum [kg m/s]",En.LM+1},
  {"LMz", "linear momentum [kg m/s]",En.LM+2},
  {"TLM", "cluster translational temperature (from CoM) [K]",&En.TLM},
  {"AMx", "angular momentum [kg m2/s]",En.AM},
  {"AMy", "angular momentum [kg m2/s]",En.AM+1},
  {"AMz", "angular momentum [kg m2/s]",En.AM+2},
  {"TAM", "cluster rotational temperature [K]",&En.TAM},
#ifdef POLAR
  {"Tpol", "POLAR DEPRECATED T of Car-Parrinello-like dipoles [K]",&En.Tpol}, /* new in 2.0c: to be used for mech.dip. */
  {"pstd","POLAR 1st iteration stderr [p.u.]",&En.polstderr},
  {"pmax","POLAR 1st iteration max. err. [p.u.]",&En.polmaxerr},
  {"plst","POLAR last iter stderr [p.u.]",&En.pollaststderr},
  {"plmx","POLAR last iter max err. [p.u.]",&En.pollastmaxerr},
  {"Pstd","POLAR stderr (vs. iterated solution, scf.test=1) [p.u.]",&En.Polstderr},
  {"Pmax","POLAR max err. (vs. iterated solution, scf.test=1) [p.u.]",&En.Polmaxerr},
  {"rate","POLAR averaged convergence rate (scf.test=1 or dV)",&scf.rate},
  {"iter","POLAR number of iterations (per MD step)",&En.polnit},
  {"X_hf","high-frequency susceptibility by direct response to ext.field",&En.chi_hf},
#endif /*# POLAR */
#if PRESSURETENSOR&PT_VIR
/* virial pressure tensor */
  {"Pvxx","site-based virial+constraint part of pressure tensor [Pa]",En.Pvir},
  {"Pvyy","site-based virial+constraint part of pressure tensor [Pa]",En.Pvir+1},
  {"Pvzz","site-based virial+constraint part of pressure tensor [Pa]",En.Pvir+2},
  {"Hz",  "Enthalpy = Epot + Ekin + p_zz V, without cutoff corrections [J/mol]",&En.Hz},
#  if PRESSURETENSOR&PT_OFF
  {"Pvyz","site-based virial+constraint part of pressure tensor [Pa]",En.Pvir+3},
  {"Pvzx","site-based virial+constraint part of pressure tensor [Pa]",En.Pvir+4},
  {"Pvxy","site-based virial+constraint part of pressure tensor [Pa]",En.Pvir+5},
#  endif /*# PRESSURETENSOR&PT_OFF */
#endif /*# PRESSURETENSOR&PT_VIR */
#if PRESSURETENSOR&PT_KIN
/* site-based kinetic pressure tensor */
  {"Pkxx","site-based kinetic part of pressure tensor [Pa]",En.Pkin},
  {"Pkyy","site-based kinetic part of pressure tensor [Pa]",En.Pkin+1},
  {"Pkzz","site-based kinetic part of pressure tensor [Pa]",En.Pkin+2},
#  if PRESSURETENSOR&PT_OFF
  {"Pkyz","site-based kinetic part of pressure tensor [Pa]",En.Pkin+3},
  {"Pkzx","site-based kinetic part of pressure tensor [Pa]",En.Pkin+4},
  {"Pkxy","site-based kinetic part of pressure tensor [Pa]",En.Pkin+5},
#  endif /*# PRESSURETENSOR&PT_OFF */
#endif /*# PRESSURETENSOR&PT_KIN */
#if PRESSURETENSOR&PT_MOL
/* CM-based kinetic pressure tensor */
  {"PKxx","molecular (CM-based) kinetic part of pressure tensor [Pa]",En.PKin},
  {"PKyy","molecular (CM-based) kinetic part of pressure tensor [Pa]",En.PKin+1},
  {"PKzz","molecular (CM-based) kinetic part of pressure tensor [Pa]",En.PKin+2},
#  if PRESSURETENSOR&PT_OFF
  {"PKyz","molecular (CM-based) kinetic part of pressure tensor [Pa]",En.PKin+3},
  {"PKzx","molecular (CM-based) kinetic part of pressure tensor [Pa]",En.PKin+4},
  {"PKxy","molecular (CM-based) kinetic part of pressure tensor [Pa]",En.PKin+5},
#  endif /*# PRESSURETENSOR&PT_OFF */
#endif /*# PRESSURETENSOR&PT_MOL */
#if (PRESSURETENSOR&PT_ANY) == PT_ANY
/* total pressure tensor */
  {"Ptxx","pressure tensor w/o cutoff corrections [Pa]",En.Ptens},
  {"Ptyy","pressure tensor w/o cutoff corrections [Pa]",En.Ptens+1},
  {"Ptzz","pressure tensor w/o cutoff corrections [Pa]",En.Ptens+2},
  {"Ptrc","P from pressure tensor w. cutoff corr: tr(Pt)/3+corr [Pa]",&En.Ptr.c},
  {"Ptrn","P from pressure tensor w/o cutoff corr : tr(Pt)/3 [Pa]",&En.Ptr.n},
#  if PRESSURETENSOR&PT_OFF
  {"Ptyz","pressure tensor w/o LJ cutoff corrections [Pa]",En.Ptens+3},
  {"Ptzx","pressure tensor w/o LJ cutoff corrections [Pa]",En.Ptens+4},
  {"Ptxy","pressure tensor w/o LJ cutoff corrections [Pa]",En.Ptens+5},
#  endif /*# PRESSURETENSOR&PT_OFF */
#endif /*# (PRESSURETENSOR&PT_ANY) == PT_ANY */
  {"Lx","x box size [p.u.=AA]",box.L},
  {"Ly","y box size [p.u.=AA]",box.L+1},
  {"Lz","z box size [p.u.=AA]",box.L+2},
  {"Jx","current density in the x direction [A/m2]",En.J},
  {"Jy","current density in the y direction [A/m2]",En.J+1},
  {"Jz","current density in the z direction [A/m2]",En.J+2},
  {"Jx0","current density in the x direction, group 0 [A/m2]",En.Jg[0]},
  {"Jy0","current density in the y direction, group 0 [A/m2]",En.Jg[0]+1},
  {"Jz0","current density in the z direction, group 0 [A/m2]",En.Jg[0]+2},
  {"Jx1","current density in the x direction, group 1 [A/m2]",En.Jg[1]},
  {"Jy1","current density in the y direction, group 1 [A/m2]",En.Jg[1]+1},
  {"Jz1","current density in the z direction, group 1 [A/m2]",En.Jg[1]+2},
  {"Jx2","current density in the x direction, group 2 [A/m2]",En.Jg[2]},
  {"Jy2","current density in the y direction, group 2 [A/m2]",En.Jg[2]+1},
  {"Jz2","current density in the z direction, group 2 [A/m2]",En.Jg[2]+2},
  {"Mx","dipole moment of the box in x (ions to centroid) [p.u.]",En.M},
  {"My","dipole moment of the box in y (ions to centroid) [p.u.]",En.M+1},
  {"Mz","dipole moment of the box in z (ions to centroid) [p.u.]",En.M+2},
  {"Mx0","raw dipole moment of the box in x, group 0 [p.u.]",En.Mg[0]},
  {"My0","raw dipole moment of the box in y, group 0 [p.u.]",En.Mg[0]+1},
  {"Mz0","raw dipole moment of the box in z, group 0 [p.u.]",En.Mg[0]+2},
  {"Mx1","raw dipole moment of the box in x, group 1 [p.u.]",En.Mg[1]},
  {"My1","raw dipole moment of the box in y, group 1 [p.u.]",En.Mg[1]+1},
  {"Mz1","raw dipole moment of the box in z, group 1 [p.u.]",En.Mg[1]+2},
  {"Mx2","raw dipole moment of the box in x, group 2 [p.u.]",En.Mg[2]},
  {"My2","raw dipole moment of the box in y, group 2 [p.u.]",En.Mg[2]+1},
  {"Mz2","raw dipole moment of the box in z, group 2 [p.u.]",En.Mg[2]+2},
  {"NHxi","Nose-Hoover degree of freedom xi = ln(s)",&En.logs},
  {"NHdx","time derivative of the Nose-Hoover degree of freedom xi = ln(s)",&En.dlogs},
  {"MTlx","MTK lambda_x scaling",En.lambda},
  {"MTly","MTK lambda_y scaling",En.lambda+1},
  {"MTlz","MTK lambda_z scaling",En.lambda+2},
  {"MTdx","time derivative of MTK lambda_x scaling",En.dlambda},
  {"MTdy","time derivative of MTK lambda_y scaling",En.dlambda+1},
  {"MTdz","time derivative of MTK lambda_z scaling",En.dlambda+2},
#ifdef RGYR
  {"Rgxx","averaged diagonalized moment of inertia of species 0",&En.RGm[0]},
  {"Rgyy","averaged diagonalized moment of inertia of species 0",&En.RGm[1]},
  {"Rgzz","averaged diagonalized moment of inertia of species 0",&En.RGm[2]},
#endif /*# RGYR */
#ifdef SLAB
  /* former xsh etc. */
  {"cx","x-drop center as determined by \"autocenter\" [AA]",&box.center[0]},
  {"cy","y-drop,trickle center as determined by \"autocenter\" [AA]",&box.center[1]},
  {"cz","z-drop,trickle,slab center as determined by \"autocenter\" [AA]",&box.center[2]},
  /* xsym etc. */
  {"symx","x-center of species <= slab.sym wrt \"autocenter\" [AA]",&En.sym[0]},
  {"symy","y-center of species <= slab.sym wrt \"autocenter\" [AA]",&En.sym[1]},
  {"symz","z-center of species <= slab.sym wrt \"autocenter\" [AA]",&En.sym[2]},
  {"Pt","tangential pressure (Ptxx+Ptyy-2*Ptzz)/3 [Pa]",&En.Pt},
  {"Ptv","virial part of tangential pressure (Pvxx+Pvyy-2*Pvzz)/3 [Pa]",&En.Ptv},
  {"Pw0","pressure from averaged forces on wall 0 [Pa]",En.Pwall},
  {"Pw1","pressure from averaged forces on wall 1 [Pa]",En.Pwall+1},
#  if SLAB & 2
  {"Wcl","dG/dsigma (rev. work) per unit interface area [Pa]",&En.Wcleave},
  {"clz0","cleave.z[0] = center of the cleaving wall 0 [Lz]",cleave.z},
  {"clz1","cleave.z[1] = center of the cleaving wall 1 [Lz]",cleave.z+1},
  {"clf0","force to cleaving wall 0 [K/AA]",En.fcleave},
  {"clf1","force to cleaving wall 1 [K/AA]",En.fcleave+1},
#  endif /*# SLAB & 2 */
#endif /*# SLAB */
#ifdef XSECTION
  {"XSx"  ,"cross-section in the x-axis",En.xsec},
  {"XSy"  ,"cross-section in the x-axis",En.xsec+1},
  {"XSz"  ,"cross-section in the x-axis",En.xsec+2},
#endif /*# XSECTION */
#ifdef COULOMB
  {"Ex"  ,"electrostatic field in x-direction [V/m]",Eext.E},
  {"Ey"  ,"electrostatic field in y-direction [V/m]",Eext.E+1},
  {"Ez"  ,"electrostatic field in z-direction [V/m]",Eext.E+2},
  {"einf","el.epsinf, dielectric constant at infinity",&el.epsinf},
  {"xinf","1/(2*el.epsinf+1), ~ constant at M^2",&el.xinf},
#endif /*# COULOMB */
  {"","",NULL}
};

void initCP(int no,int nbit,char *col4, char *col5,int corr) /******* initCP */
/*
   nbit=0 : normal write in float[NCP] format
   nbit>0 : packed, nbit bits
*/
{
  float *CPheader=NULL;
  char key;
  FILE *CP;
  int ncp=0,lengthchanged=0;
  char *corrinfo = corr&1 ? "incl. LJ cutoff corr. "
                          : "w/o LJ cutoff corr. ";
#ifdef ANCHOR
  int anchorincpi=0;
#endif /*# ANCHOR */

  if (!option('e')) return;

  if (!CPinfo) {
    /* this initialization (of CPinfo and ssd0->..) is done only once */
    FILE *CPI=fopen(Fn("cpi"),"rt");
    char line[LINELEN];
    int nssd=0,i,col,e,CPinfo0len;
    struct ssd_s *ssd,**next=&ssd0;
    char *c,*info="?";

    if (nbit<0 || abs(nbit-3)<3 || nbit>24) {
      ERROR(("CPnbit=%d is outside range 6..24",nbit))
      nbit=0; }
    CPnbit=nbit;

    if (sizeof(CPinfo0)%4!=1) ERROR(("\"%s\" is bad CPinfo0",CPinfo0))

    strcpy(line,col5); strcat(line,"   "); memcpy(CPinfo0+8,line,4);
#ifdef LOG
    if (No.first) {
      char x[8];
      sprintf(x,"%d  ",No.first);
      copy(CPinfo0+5,x,3); }
#endif /*# LOG */
    strcpy(line,col4); strcat(line,"   "); memcpy(CPinfo0+4,line,4);

    /* whether current density or dipole moment monitored in column 8 */
    key=No.ion?'J':'M';
    if (Eext.isE<0) copy(CPinfo0+20,string("%c%c  ",key,Eext.isE+3+'x'),4);
    else copy(CPinfo0+20,string("|%c| ",key),4);

#if 0
    /*
       THIS IS WRONG because does not match the order of recordCP in maincps.c
       Column 8 will be active for COULOMB, but it may be zero
       unless lag.J/lag.M is specified.
       Planned remedy: make M,J dynamic (in .cpi) ?
    */
#  if COULOMB>=0 /* Fennel-Gezelter ? */
    if (lag.J==0 && lag.M==0) CPinfo0[20]=0;
#  endif /*# COULOMB>=0 */
    // #endif /*# COULOMB */
#  ifdef FREEBC
    if (lag.J==0 && lag.M==0) CPinfo0[20]=0;
#  endif /*# FREEBC */
#endif /* wrong */ /*# 0 */

    underline("convergence profile");
    prt("HINT: use \"showcp [OPTIONS] %s.cp\" to show and analyze",simils.simname);
    prt("INFO: if no unit given, p.u.=program unit assumed, see section \"units\" above");
    prt(" 1:Etot total energy or extended Hamiltonian [K]");
    prt(" 2:T    kinetic temperature (cf. Tin, Ttr below) [K]");
    col=2;

    //    fprintf(stderr,"CPinfo0=\"%s\" len=%d\n",CPinfo0,(int)strlen(CPinfo0));

    /* keywords in CPinfo0 (starts with Epot, composed via #if) -- cumbersome! */
    for (c=CPinfo0; *c>' '; c+=4) {
      col++;
      prt_("%2d:%c%c%c%c ",col,c[0],c[1],c[2],c[3]);
      if (!memcmp("Epot",c,4)) prt("total potential energy %s[J/mol]",corrinfo);
      else if (!memcmp("Tin ",c,4)) prt("temperature from internal motions and rotations [K]");
      else if (!memcmp("Ttr ",c,4)) prt("translational temperature from molecular centers of mass [K]");
// prior V2.8c: else if (!memcmp("Ein ",c,4)) prt("intramolecular potential energy [J/mol]");
      else if (!memcmp("Ep+k",c,4)) prt("potential+kinetic energy (w/o Nose) [J/mol]");
      /* this column changed via initCP */
      else if (!memcmp("rho ",c,4)) prt("density [kg/m3]");
      else if (!memcmp("svdW",c,4)) prt("van der Waals sigma of the i-j pair being adjusted [AA]");
      else if (!memcmp("PdVm",c,4)) prt("pressure by CM-based virtual volume change %s[Pa]",corrinfo);
      else if (!memcmp("PdVa",c,4)) prt("pressure by atom-based virtual volume change %s[Pa]",corrinfo);
      else if (!memcmp("P   ",c,4)) prt("pressure (from virial and kin. energy) %s[Pa]",corrinfo);
      else if (!memcmp("Pref",c,4)) prt("pressure (as selected by variable virial) %s[Pa]",corrinfo);
      else if (!memcmp("|M| ",c,4)) prt("sqrt(<M^2>), simulation cell dipole moment [p.u.]; cf. lag.M");
      else if (!memcmp("Mx  ",c,4)) prt("<Mx>, simulation cell dipole moment in x [p.u.]");
      else if (!memcmp("My  ",c,4)) prt("<My>, simulation cell dipole moment in y [p.u.]");
      else if (!memcmp("Mz  ",c,4)) prt("<Mz>, simulation cell dipole moment in z [p.u.]");
      else if (!memcmp("|J| ",c,4)) prt("sqrt(<J^2>), simulation cell current density [A/m2]; cf. lag.J");
      else if (!memcmp("Jx  ",c,4)) prt("<Jx>, simulation cell current density in x [A/m2]");
      else if (!memcmp("Jy  ",c,4)) prt("<Jy>, simulation cell current density in x [A/m2]");
      else if (!memcmp("Jz  ",c,4)) prt("<Jz>, simulation cell current density in x [A/m2]");
#ifdef SHEAR
      else if (!memcmp("Cvel",c,4)) prt("shear-induced velocity [p.u.]");
#endif /*# SHEAR */
#ifdef POLAR
      else if (!memcmp("self",c,4)) prt("self-field [p.u.]");
#endif /*# POLAR */
#ifdef RGYR
      else if (!memcmp("Rgyr",c,4)) prt("radius of gyration [AA]");
      else if (!memcmp("ende",c,4)) prt("end-to-end distance [AA], see rg.end[]");
#endif /*# RGYR */
#ifdef XSECTION
      else if (!memcmp("Xsec",c,4)) prt("cross section [AA2]");
#endif /*# XSECTION */
      /* to be extended for more versions -- or rehacked ... */
      else _n }

    i=0;
    if (CPI) {
      while (fgets(line,LINELEN,CPI)) {
        i++;
        if (line[0]!='!' && line[0]!='#' && line[0]!='\n') {
          alloczero(ssd,sizeof(*ssd));
          ssd->stat=line[0]=='+';
          e=sscanf(line+ssd->stat,"%d%d%s",&ssd->indx[0],&ssd->indx[1],ssd->name);

          if (e>=2) {
            /* recording atom-atom distance */
            if (e==2) sprintf(ssd->name,"ss%d",nssd+1);
            loop (e,0,2) if (ssd->indx[e]<0 || ssd->indx[e]>=No.s) {
              ERROR(("%s, line %d: atom index %d out of range [0,%d)",
                     lastFn,i,ssd->indx[e],No.s))
              ssd->indx[e]=0; }
            prt("%2d:%-4s |%d-%d| distance [AA]",
              col+i,ssd->name,ssd->indx[0],ssd->indx[1]); }

          else if (e==0) {
            /* check names from the list CPItab */
            e=sscanf(line+ssd->stat,"%s",ssd->name);
            ssd->indx[0]=-9; /* no unit conversion */

            if (e) {
              int i;
#ifdef ANCHOR
              if (option('k')<0 && (-option('k'))&2) loop (i,0,anchor.col) {

                if (!strcmp(anchor.rec[i].info,ssd->name)) {
                  anchorincpi++;
                  info="see ANCHOR block above";
                  ssd->u.q=&anchor.rec[i].var;
                  goto CPIfound; } }
#endif /*# ANCHOR */

              for (i=0; CPItab[i].var; i++) {

                if (!strcmp(CPItab[i].name,ssd->name)) {
                  info=CPItab[i].info;
                  ssd->u.q=CPItab[i].var;

                  /* NEW: control changed, measurement duplicated, variable measuredrift removed
                     the following names in .cpi removed: |CM| CMdr |lM| lMdr |aM| aMdr */

                  /* denote conversion of p.u. to J/mol, V/m, Pa on .cp */
                  if (strstr(CPItab[i].info,"[J/mol]")) ssd->indx[0]=-1;
                  if (strstr(CPItab[i].info,"[V/m]")) ssd->indx[0]=-2;
                  if (strstr(CPItab[i].info,"[Pa]")) ssd->indx[0]=-3;

                  goto CPIfound; } }
              e=0;
              ERROR(("%s: keyword \"%s\" is not known, see struct CPItab in simmeas.c",lastFn,ssd->name))
              CPIfound:; }
            prt("%2d:%-4s %s", col+i, ssd->name, info); }
          else {
            e=0;
            ERROR(("%s: Bad format or not enough data in line %d.",lastFn,i)) }

          if (e) {
            nssd++;
// padding removed!!!
//            while (strlen(ssd->name)<4) strcat(ssd->name," ");
            ssd->next=NULL;
            *next=ssd;
            next=&ssd->next; }
        } /* noncomment line */
      } /* while line of CPI */

      fclose(CPI);

#ifdef ANCHOR
      if (option('k')<0 && (-option('k'))&2)
        if (anchorincpi==0) WARNING(("ANCHOR: -k-2 selected but no keyword found in %s",lastFn))
#endif /*# ANCHOR */
      prt("%s: %d additional variables and site-site distances recorded",lastFn,nssd);
      CPinfo0len=strlen(CPinfo0);
      if (CPinfo0len&3) ERROR(("internal: CPinfo0 not multiple of 4 bytes"))
      alloczero(CPinfo,CPinfo0len+4*nssd+1); // !V3.6l +1 because trailing 0 needed since NCP is derived from length
      copy(CPinfo,CPinfo0,CPinfo0len);

      for (ssd=ssd0,c=CPinfo+CPinfo0len; ssd; ssd=ssd->next,c+=4) {
        /* keyword padded by SPACE here */
        memset(c,' ',4);
        if (strlen(ssd->name)>4) {
          copy(c,ssd->name,3);
          c[3]=strlast(ssd->name)[0]; }
        else
          copy(c,ssd->name,strlen(ssd->name)); } }
    else {
      CPinfo=CPinfo0;
      prt("File %s is missing, only default items will be recorded.",lastFn); } }

  NCP=2+strlen(CPinfo)/4;

  if (option('e')>0) {
    if (option('e')==1) {
      option('e')=2;
      WARNING(("option -e1 not allowed, changed to -e2")) }
    Min(NCP,option('e'))
    prt(">>> because of option -e, only %d columns recorded",option('e')); }

  if (init<2) {

    if (!(CP=fopen(Fn(CPnbit?"cpz":"cp"),"rb")))
      ERROR(("Cannot open %s (CPnbit=%d) for reaading",Fn(CPnbit?"cpz":"cp"),CPnbit))
    else {
      prt("Convergence profile file %s opened for reading",lastFn);
      CPheader=getcpheader(CP,&ncp,option('v')&0xfff8);

      if (ncp!=NCP) {
        lengthchanged++;
        WARNING(("%s contains %d fields but %d requested: extra fields %s",
                 lastFn,ncp,NCP,ncp<NCP?"truncated":"padded by 0"))
        if (ncp>NCP) {
          float *newCPheader;
          allocarray(newCPheader,ncp);
          memset(newCPheader+NCP,'-',(ncp-NCP)*4);
          memcpy(newCPheader,CPheader,8);
          free(CPheader);
          CPheader=newCPheader; }
        NCP=ncp; }
      fclose(CP); } }

  prt("CP record %d floats:",NCP);

  if (no && NCP>=2) {
    if ((long)no*sizeof(float)*NCP>(long)AllocSizeLim)
      ERROR(("Convergence profile data (%g GiB) longer than limit %g GiB,\n\
*** decrease number of cycles (no) or increase AllocSizeLim",
             (long)no*sizeof(float)*NCP/(double)0x40000000,
             AllocSizeLim/(double)0x40000000))
    ralloc(CPbuf,(int4)no*sizeof(float)*NCP); }

  prtCPinfo(CPinfo);
  if (CPheader) {
    if (!lengthchanged) if (memcmp(CPinfo,CPheader+2,NCP*4-8)) {
        prt("\
raw header=\"%s\"\n\
   in file=\"%s\"",CPinfo,CPheader+2);
        prt("header in file %s:",lastFn);
        prtCPinfo((char*)(CPheader+2));
        WARNING(("inconsistent header of CP file\n\
*** this is most often caused by using init=\"cont\" (0) instead of \"start\" (2)")) }
    free(CPheader); }

  prt("NB: for showcp -a2/-p2, column 1 = time and other columns +1\n\
CP conversion factors:\n\
  column 0 (line number, with showcp -b1) to s:           %.9g*n\n\
  column 1 (total energy in K) to J/mol:                  %.9g*A\n\
  column 2 (temperature in K) to kinetic energy in J/mol: %.9g*B\n\
  column 4 (density in kg/m3) to volume in AA3:           %.9g/D\n",
      No.dt*timeunit,Eunit,No.f*(Eunit/2),rhounit*No.mass);
}

void saveCP(int icyc,time_t stoptime) /****************************** saveCP */
{
  static char CPmode[3]="wb";
  static FILE *CP;
  float hdr=CPmark;
  char *aux;
  int4 ncp;

  if (!CPbuf) return;

  if (sizeof(int4)!=sizeof(float)) ERROR(("sizeof(int4)!=sizeof(float))"))
  if (sizeof(float)!=4) ERROR(("sizeof(float)!=4"))

  aux=CPnbit?"cpz":"cp";
  if (init<=1)
    CPmode[0]='a';
  else {
    CPmode[0]='w';
    backup(aux); }

  CP=fopen(Fn(aux),CPmode);
  if (!CP) ERROR(("%s: open",lastFn))
  if (CPmode[0]=='a')
    prt("appending %s at position %ld", lastFn, ftell(CP));
  if (CPmode[0]=='w') {
    /* create header */
    fwrite(&hdr,sizeof(float),1,CP);
    ncp=NCP;
    fwrite(&ncp,sizeof(int4),1,CP);
    fwrite(CPinfo,4,NCP-2,CP);  }
  /* here we are after timestamp or header */
  /* new: info on time and timestep */
  if (NCP>=5) {

    double rec[2]={No.t0,No.dt}; /* to be compatible with float cook */

    hdr=CPmarkW;
    fwrite(&hdr,sizeof(float),1,CP);
    fwrite(rec,sizeof(double),2,CP);
    fwrite(CPinfo+3,4,NCP-5,CP);  }
  else if (NCP>=3) {
    double rec=t; /* well, to be compatible with float */

    hdr=CPmarkU;
    fwrite(&hdr,sizeof(float),1,CP);
    fwrite(&rec,sizeof(double),1,CP);
    fwrite(CPinfo+1,4,NCP-3,CP);

    rec=No.dt;
    hdr=CPmarkU;
    fwrite(&hdr,sizeof(float),1,CP);
    fwrite(&rec,sizeof(double),1,CP);
    fwrite(CPinfo+1,4,NCP-3,CP); }

  if (CPnbit) appendpakcp(CP,NCP,CPnbit,(void*)CPbuf,icyc);
  else fwrite(CPbuf,sizeof(float)*NCP,icyc,CP);

  /* ending record with time stamp */
  memset(CPbuf,0,NCP*sizeof(float));
  CPbuf[0]=CPmarkT;

  /* changed in V3.3m, 11/2019 */
  if (NCP>=3) {
    /* sizeof(time_t)=8 assumed */
    if (NCP>=5) /* very old legacy position */
      copy(CPbuf+NCP-sizeof(stoptime)/sizeof(float),&stoptime,sizeof(stoptime));
    copy(CPbuf+1,&stoptime,sizeof(stoptime)); }
  else if (NCP==2) {
    unsigned u=stoptime;
    copy(CPbuf+1,&u,sizeof(float)); }

#if 0
  /* "new" version abandoned 2019 */
  aux=ctime((void*)&stoptime);
  l=min(strlen(aux),sizeof(float)*(NCP-1));
  copy(CPbuf+1,aux,l);
  aux=(char*)(CPbuf+1);
  while (l--) { if (*aux=='\n') *aux=0; aux++; }
  ((char*)(CPbuf+NCP))[-1]=0; /* ??? */
#endif /*# 0 */

  fwrite(CPbuf,sizeof(float)*NCP,1,CP);

  prt("closing %s at position %ld", lastFn,ftell(CP));
  if (fclose(CP)) if (!CP) ERROR(("%s: close",lastFn))
}


/*
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  site-site radial distribution functions
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

static int global_rdfgrid;

void initrdf(double gr,double cutoff,char *fn) /******************** initrdf */
/*
  gr=grid (subintervals/length unit)
  initializing histograms for rdf's,
  gr=0: none
  gr<0: as |gr| but all site-site histograms merged
  assigned to ss->rdf in setss which is called by initss
*/
{
  int i,j,nss=0;
  int sp;
  int *sitecount;
  int nh=(int)(cutoff*abs(gr)+1.0000001);
  int size=sizeof(rdf_t)-2*sizeof(int4)+sizeof(int4)*nh;
  rdf_t *r;
  FILE *ff=NULL;

  struct ss_s {
    struct ss_s *next;
    char s1[8],s2[8];
  } *ss,*ss0=NULL;

  underline("nonbonded pair terms and radial distribution functions");
  if (!gr) prt("Radial distribution functions are turned off.");

  if (cutoff>8.9e9) cutoff=100; /* FREEBC, NIBC */

  if (cutoff<0 || cutoff>1e4) ERROR(("rdf.cutoff=%g out of range [0..1e4]",cutoff))
  if (nh>100000) ERROR(("number of RDF histogram bins=%d is too large",nh))

  if (gr>0 && !(ff=fopen(fn,"rt")))
    prt("%s not found, RDFs will be recorded for all pairs of site types.",fn);

  if (ff) {
    char line[LINELEN],*s1,*s2;

    prt("reading %s: site type pairs",fn);

    while (fgets(line,LINELEN,ff)) if (!strchr("!#",line[0])) {
      s1=strtok(line," \t\n"); if (!s1) continue;
      if (s1[0]=='*') {
        if (nss) ERROR(("%s: wildcard specifying all pairs comes after %d pairs.",lastFn,nss))
        nss=0; ss0=NULL; /* memory leak */
        fclose(ff); ff=NULL;
        goto escape; }
      if (strlen(s1)==1 && (islower(s1[0])||isdigit(s1[0]))) continue;
      s2=strtok(NULL," \t\n");
      if (!s2) { ERROR(("%s: %s what?",fn,s2)) continue; }

      alloczero(ss,sizeof(struct ss_s));
      copy(ss->s1,s1,5);
      copy(ss->s2,s2,5);
      ss->next=ss0; ss0=ss;
      nss++; }
    prt("%d pairs for rdf selected",nss); }

 escape:

  global_rdfgrid=gr;

  /* now site-type statistics (option -j ignored) */

  allocarrayzero(sitecount,nsites);

  loop (sp,0,nspec) {
    specinfo_t *s=spec[sp];
    loop (i,0,s->ns) sitecount[s->si[i].st] += s->N; }

  ralloc(rdf,nsites*sizeof(rdf_t**));

  loop (i,0,nsites) {
    ralloc(rdf[i],nsites*sizeof(rdf_t*));
    loop (j,0,nsites) rdf[i][j]=NULL; }

  if (option('v')&1) {
    prt("\n _i_site_count  _i_site_count  _i_site_count  _i_site_count  _i_site_count");
    loop (i,0,nsites)
      prt_("%3d %-4s %4d%c",i,sitedef[i].name,sitecount[i],"    \n"[i%5]);
    _n }

  if (gr) loop (i,0,nsites) {

    loopto (j,0,i) {
      if (gr>0 || (i|j)==0) {
        if (ff) {
          for (ss=ss0; ss; ss=ss->next) if (
             (!strcmp(ss->s1,sitedef[i].name) && !strcmp(ss->s2,sitedef[j].name))
           ||(!strcmp(ss->s1,sitedef[j].name) && !strcmp(ss->s2,sitedef[i].name)) )
            goto incl;
          r=NULL;
          goto noincl; /* not include this pair */ }
      incl:
        sdsralloczero(r,size); /* cleared all .. purify! */
        r->nmeas=r->npair=0;
        r->indx[0]=i; r->indx[1]=j;
        if (gr<=0) {
          r->ns[0]=r->ns[1]=No.s;
          strcpy(r->name[0],"*");
          strcpy(r->name[1],"*"); }
        else {
          r->ns[0]=sitecount[i];
          copy(r->name[0],sitedef[i].name,
               min(sizeof(r->name[0]),sizeof(sitedef[0].name)));
          r->ns[1]=sitecount[j];
          copy(r->name[1],sitedef[j].name,
               min(sizeof(r->name[0]),sizeof(sitedef[0].name))); }
        r->grid=(double)abs(gr);
        r->V=0;
        r->nhist=nh; }
      else
        r=rdf[0][0];
    noincl:
      rdf[i][j]=rdf[j][i]=r; } }

  free(sitecount);
  if (ff) fclose(ff);
  for (ss=ss0; ss; ) {
    struct ss_s *new=ss->next;
    free(ss);
    ss=new; }
}

void advancerdf(void) /****************************************** advancerdf */
{
  int i,j;
  rdf_t *r;

  if (global_rdfgrid && measure)
    loop (i,0,nsites) loopto (j,0,i) {
      r=rdf[i][j];
      if (r) { r->nmeas++; r->V+=box.L[0]*box.L[1]*box.L[2]; }
      if (global_rdfgrid<0) return; }
}

void loadrdf(double rdfgrid) /************************************** loadrdf */
{
  int ng = rdfgrid<0 ? 1 : nsites;
  int nsi,i,j;
  double gr;

  if (!rdfgrid) return;

  VarOpen(Fn("rdf"),"r");
  VarRead(&gr,sizeof(gr));
  VarRead(&nsi,sizeof(nsi));
  if (gr!=rdfgrid || nsi!=nsites)
    ERROR(("%s: grid or # of sites inconsistent",lastFn))
  loop (i,0,ng)
    loopto (j,0,i)
      if (rdf[i][j]) {
        VarReadSds(rdf[i][j]);
        if (rdf[i][j]->nhist*sizeof(rdf[i][j]->hist[0])+((char*)rdf[i][j]->hist-(char*)rdf[i][j]) != rdf[i][j]->size)
        ERROR(("Inconsistent RDF file. If you continue with RDF obtained by cook < V2.7f,\n\
*** use  rdfconv %s  first",lastFn)) }

  VarClose();
}

void saverdf(double rdfgrid) /************************************** saverdf */
{
  int ng = rdfgrid<0 ? 1 : nsites, i,j;

  if (!rdfgrid) return;

  backup("rdf");
  VarOpen(Fn("rdf"),"w");
  VarPut(&rdfgrid,sizeof(rdfgrid));
  VarPut(&nsites,sizeof(nsites));
  loop (i,0,ng) loopto (j,0,i)
    if (rdf[i][j]) VarPutSds(rdf[i][j]);
  VarClose();
}

#if PARALLEL==1
void parallelizerdf(int replicate)
/*
  replicate=1: replicate the histograms to be able to sum them concurrently
  replicate=0: sum back to one histogram
*/
{
  int ih,ith,i,j;

  loop (i,0,nsites) loopto (j,0,i) if (replicate) {
    /* replicate histograms, grid, nhist (more not needed) */
    rallocarrayzero(sstab[i][j].rdfs,No.th);
    sstab[j][i].rdfs=sstab[i][j].rdfs;
    loop (ith,0,No.th)
      if (!rdf) sstab[i][j].rdfs[ith]=NULL; /* should not happen */
      else if (No.measureserial)
        sstab[i][j].rdfs[ith]=sstab[i][j].rdf; /* may be NULL */
      else if (sstab[i][j].rdf) {
        if (ith) {
          /* replica */
          sdsralloczero(sstab[i][j].rdfs[ith],sstab[i][j].rdf->size);
          sstab[i][j].rdfs[ith]->grid=sstab[i][j].rdf->grid;
          sstab[i][j].rdfs[ith]->nhist=sstab[i][j].rdf->nhist; }
        else {
          /* the original structure */
          sstab[i][j].rdfs[ith]=sstab[i][j].rdf; } } }
  else {
    /* sum back */
    if (sstab[i][j].rdf && !No.measureserial)
      loop (ih,0,sstab[i][j].rdf->nhist)
        loop (ith,1,No.th)
          sstab[i][j].rdf->hist[ih]+=sstab[i][j].rdfs[ith]->hist[ih];
    }
}
#endif /*# PARALLEL==1 */

/*
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                           dihedral angle distribution
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifdef DIHHIST

static int diheq(char type1[4][6],char type2[4][6])
/*
  returns 1 if the dihedral has already been included
  (according to the rules given by #define DIHHIST)
*/
{
#  define EQ(I,J) !strcmp(type1[I],type2[J])

#  if DIHHIST==2
  return (EQ(1,1) && EQ(2,2)) || (EQ(1,2) && EQ(2,1));
#  elif DIHHIST==4
  return
      (EQ(0,0) && EQ(1,1) && EQ(2,2) && EQ(3,3))
   || (EQ(0,3) && EQ(1,2) && EQ(2,1) && EQ(3,0));
#  elif DIHHIST==-1
  return 0; /* all included separately */
#  elif DIHHIST==1
  return 1; /* all dihedrals summed */
#  else /*#!DIHHIST==2!DIHHIST==4!DIHHIST==-1!DIHHIST==1 */
#    error "illegal DIHHIST"
#  endif /*#!DIHHIST==2!DIHHIST==4!DIHHIST==-1!DIHHIST==1 */
#  undef EQ
}

dihhist_t *dihhead;

void includedihhist(torsion_t *t, char type[4][6], int sp) /* includedihhist */
/*
 includes a dihedral provided that it has not been included yet (the rules
 when a dihedral is considered as the same depend on #define DIHHIST)
*/
{
  dihhist_t *dih;

  for (dih=dihhead; dih; dih=dih->next)
    if (diheq(type,dih->type)) goto already;
  alloczero(dih,sizeof(dihhist_t));
  dih->next=dihhead;
#  if DIHHIST==-1
  dih->sp=sp;
#  else /*# DIHHIST==-1 */
  dih->sp=-1; /* species irrelevant */
#  endif /*#!DIHHIST==-1 */
  dihhead=dih;
  copy(dih->type,type,sizeof(dih->type));
 already:
  t->parm.dihhist=dih;
  if (abs(t->parm.n)==3) dih->gauche_trans=2*PI/3;
  else dih->gauche_trans=PI/2;
}

void initdihhist(int dihgrid) /********************************* initdihhist */
{
  dihhist_t *dih;
  double q=dihgrid/(2*PI+No.eps);

  for (dih=dihhead; dih; dih=dih->next) {
    dih->grid=dihgrid;
    dih->q=q;
    if (dihgrid) ralloczero(dih->hist,sizeof(unsigned)*dihgrid);
    else dih->hist=NULL; }
}

double gauchetrans(int key) /*********************************** gauchetrans */
/*
  key=-3: assigns all gauche[cis]/trans statistics to 0
  key=-2: adds gauche/trans statistics to StaAdd and assigns 0
  key=-1: calculates total gauche ratio
  key>0: returns gauche/all for dihedral# key
         (NEW: key==0 not alloweed)
*/
{
  char name[20];
  dihhist_t *dih;
  int sum,idih=0;
  int sumg=0,sumd=0;

  for (dih=dihhead; dih; dih=dih->next) {

    sum=dih->trans+dih->gauche;

    if (key==-2 && sum) {
      sprintf(name,"%s(%d)",dih->gauche_trans<2?"cis":"gauche",idih+1);
      StaAdd(name,dih->gauche/(double)sum); }

    if (key<-1) dih->gauche=dih->trans=0;
    else if (key==-1) {
      sumg+=dih->gauche;
      sumd+=sum; }
    else if (idih+1==key) return dih->gauche/(double)sum;

    idih++; }

  if (key>=0) ERROR(("no dihedral # %d",key))

  return sumg/(sumd+1e-33);
}

void printdihhist(void) /************************************** printdihhist */
{
  dihhist_t *dih;
  unsigned sum,trans,gauche;
  int i,idih=0;
  FILE *f=NULL;
  double d;

  if (!dihhead) return;

  for (dih=dihhead; dih; dih=dih->next) {
    if (!dih->grid) continue;

#  if 0
    if (dih->hist[dih->grid]) {
      if (option('v')&1) prt("next dihedral hist[last]=%d",dih->hist[dih->grid]);
      dih->hist[dih->grid-1]+=dih->hist[dih->grid];
      dih->hist[dih->grid]=0; }
#  endif /*# 0 */

    trans=gauche=0;
    loop (i,0,dih->grid)
      if (i<2*dih->grid/3) gauche+=dih->hist[i];
      else trans+=dih->hist[i];
    sum=trans+gauche;

    if (sum || option('v')&2)
      prt_("dihedral (%d) %s-%s-%s-%s : ",
           idih+1,
           dih->type[0],dih->type[1],dih->type[2],dih->type[3]);

    if (sum) {
      if (!f) {
        _n
        backup("dia");
        prt_("                            ");
        f=fopen(Fn("dia"),"wt"); }

      fprintf(f,"# (%d) %s-%s-%s-%s : ",
              idih+1,
              dih->type[0],dih->type[1],dih->type[2],dih->type[3]);
      fprintf(f,"%d points, %.5f trans, %.5f %s\n",
              sum,(double)trans/sum,(double)gauche/sum,
              dih->gauche_trans<2?"cis":"gauche"); }

    if (sum || option('v')&2) prt_("%u points",sum);

    if (sum) prt(", %.5f trans, %.5f %s",
                 (double)trans/sum,(double)gauche/sum,
                 dih->gauche_trans<2?"cis":"gauche");
    else if (option('v')&2) _n

    if (sum) {
      loop (i,0,dih->grid) {
        d=dih->hist[i]/((double)sum/dih->grid*(2*PI));
        /*        if (i==0) fprintf(f,"%5.1f %7.5f\n",0.0,d); */
        fprintf(f,"%5.1f %7.5f\n", (i+0.5)*360.0/dih->grid,d); }
      /*      fprintf(f,"%5.1f %7.5f\n",360.0,d);*/
      fprintf(f,"\n"); }
    if (sum) idih++; }

  if (f) fclose(f);
}

void loaddih(int dihgrid) /***************************************** loaddih */
{
  dihhist_t *dih;
  struct {
    double dummy1;
    dihhist_t dih;
    double dummy2;} aux;

  if (!dihgrid) return;

  VarOpen(Fn("dih"),"r");
  for (dih=dihhead; dih; dih=dih->next) {
    VarRead((char*)aux.dih.type,VarFile.size); /* note: char aux.dih.type[][] */
    if (aux.dih.grid!=dih->grid
      || memcmp(aux.dih.type,dih->type,sizeof(aux.dih.type)))
      ERROR((
             "%s inconsistent (grid=%d type=%s, should be %d %s - try DIHOLDFMT?)",
             Fn("dih"),aux.dih.grid,aux.dih.type,dih->grid,dih->type))
        /* OLD version (grid in 180, safe bound): (dih->grid+1) */
        VarRead(dih->hist,dih->grid*sizeof(dih->hist[0])); }
  VarClose();
}

void savedih(int dihgrid) /***************************************** savedih */
{
  dihhist_t *dih;

  if (!dihgrid) return;

  backup("dih");
  VarOpen(Fn("dih"),"w");
  for (dih=dihhead; dih; dih=dih->next) {
    VarPut(dih->type,sizeof(dihhist_t)-(int)((char*)dih->type-(char*)dih));
    VarPut(dih->hist,dih->grid*sizeof(unsigned)); }
  VarClose();
}
#endif /*# DIHHIST */

static double chi2eps(double chi)
/* "b.c.-dependent susceptibility -> dielectric constant */
{
  return 1+1/(1/chi-1/(1+2*el.epsinf));
}

void prtdielconst(double V,int ewald) /************************ prtdielconst */
/*
  ewald=1: ewald-based
  ewald=0: meant for short-range potentials - restricted
*/
{
  double sat,factor,T;
  vector VarM={0,0,0};
  double sumM,sumMM,MM,dMM,xM=0;
#ifdef POLAR
  double chi_alpha=0;
#endif /*# POLAR */
  int sp,i,ichi,nchi,Nm;
  char *MMstring=ewald?"EwM^2":"M^2";
  char *Mfmt=ewald?"EwM%c":"M%c";
  struct chi_s {
    double M; /* from var M (low-frequency) */
    double dM; /* error */
    double hf; /* high-frequency reponse (POLAR only) */
    double dhf; /* error (POLAR only) */
    double eps_hf; /* high-frequency diel. const. (POLAR only) */
    int Nm;     /* number of measurements */
  } chi[4]; /* 0=nonpolar, or with hf from polarizability (by Clausius-Mossotti)
               1=as 0 with error from 3 components
               2=as 0 with hf measured
               3=M,dM contain both low and high frequency (0+2 merged) */

  if (ewald && Eext.isE && Eext.f) {
    int k;

    underline("oscillating external electric field");

    prt("frequency = %g [GHz]",Eext.f/timeunit);
    loop (k,0,DIM) if (Eext.E[k]) {
      char *c=string("EwM%c*cos",k+'x');
      double cosav=StaMean(c),coser=StaError;
      char *s=string("EwM%c*sin",k+'x');
      double sinav=StaMean(s),siner=StaError;

      if (StaError)
        WARNING(("Fourier components of EwM%c not found",k+'x'))
      else
        prt("E%c=%g=%g V/m  phase=%g  <%s>=%g %g  <%s>=%g %g",
            k+'x',
            Eext.E[k],
            Eext.E[k]/(chargeunit/forceunit),
            Eext.phase[k],
            c,cosav,coser,
            s,sinav,siner); }
    return; }

  underline("dielectric constant");
  /* sum of dipole moments */

  prt("Ewald summation %sused:\n\
- the dipole moments have been accumulated every %s\n\
- the calculations will be based on variables denotes %sM^2 etc.",
      ewald?"":"not ",
      ewald?"step":"cycle",
      ewald?"Ew":"");
  if (ewald)
    prt("Molecule dipole moments %savailable because lag.M=%d",lag.M?"not ":"",lag.M);

  sumM=sumMM=0;
  loop (sp,0,nspec) {
#ifdef POLAR
    sumM+=StaMean(string("mu_perm[%d]",sp))*spec[sp]->N;
    if (StaError) sumM+=StaMean(string("dipmom[%d]",sp))*spec[sp]->N; /* legacy name */
    sumMM+=StaMean(string("mu_tot[%d]",sp))*spec[sp]->N;
    if (StaError) sumMM+=StaMean(string("dipmomboth[%d]",sp))*spec[sp]->N; /* legacy name */
#else /*# POLAR */
    sumM+=StaMean(string("mu[%d]",sp))*spec[sp]->N;
    if (StaError) sumM+=StaMean(string("dipmom[%d]",sp))*spec[sp]->N; /* legacy name */
#endif /*#!POLAR */
  }
  if (!StaError) {
    prt("sum of permanent dipole moments:\n\
  M_perm = %g p.u. = %g D = %g C m",
      xM=sumM,sumM*Debye,sumM*(chargeunit*lengthunit));
#ifdef POLAR
    prt("sum of permanent+induced dipole moments\n\
  M_tot = %g p.u. = %g D = %g C m",
        xM=sumMM,sumMM*Debye,sumMM*(chargeunit*lengthunit));
#endif /*# POLAR */
  } else {
    static int pass=0;
    if (0==pass++)
      WARNING(("In order to report saturation, the dipole moments of individual\n\
*** molecules should be measured. Use at least lag.M=1 to do this.")) }

  Fn("cpa");
  prt("Hint: |M| or Mz is usually in column 8.\n\
To show saturation after showcp has generated %s, run:\n\
  plot %s:0:H/%f\n",lastFn,lastFn,xM);

  T=StaMean("Tkin");
  if (StaError) {
    WARNING(("no temperature, eps skipped"))
    return; }
  factor=4*PI/(3*T*V);

  /* this is active with COULOMB and COULOMB<0 .. what about Fennel-Gezelter? */
  MM=StaMean(MMstring);

  if (StaError) {
    WARNING(("<%s> not available, eps skipped",MMstring))
    return; }

  dMM=StaStdErr(MMstring);
  Nm=StaN(MMstring);

  prt("MM = <%s> = %g %g (in p.u.)",MMstring,MM,dMM);

  if (MM==0) {
    prt("Not dipolar system, calculation of dielectric constant skipped",MMstring);
    return; }

  prt("\nFLUCTUATION ROUTE");
  prt("These calculations are based on a fluctuation formula with el.epsinf=%g",el.epsinf);
  if (Eext.isE) prt("WARNING: EXTERNAL ELECTRIC FIELD");
  if (Eext.isB) prt("WARNING: EXTERNAL MAGNETIC FIELD");

  if (sumM) prt("Saturation: stdev(MM)/M_perm = %g",sqrt(MM)/sumM);
#ifdef POLAR
  if (sumM) prt("Saturation: stdev(MM)/M_tot = %g",sqrt(MM)/sumMM);
  if (sumMM) Eext.sumM=sumMM; else
#endif /*# POLAR */
    Eext.sumM=sumM;

  if (Eext.sumM) {
    sat=sqrt(MM)/Eext.sumM;
    if (sat>0.1) prt("*** WARNING: The dipole moment of the cell is %.3f of the saturated one.\n\
*** The dielectric constant may be affected by nonlinear response.\n\
*** Consider decreasing el.epsinf or increasing the box size.",sat); }
  else
    prt("WARNING: saturation not known, consider lag.M");

  sat=MM;
  /* subtracting <Mx>^2 etc. if needed, check of saturation */
  if (Eext.isE) {
    loop (i,0,3) if (Eext.E[i]) {
      char *Mc=string(Mfmt,i+'x');
      double av=StaMean(Mc);

      if (StaError) {
        WARNING(("With ext. field: <M%c> not found, skipped",i))
        break; }
      prt("<M%c>^2=%g^2 subtracted from <M^2> because of the elst field",i+'x',av);
      MM-=av*av; }
    prt("corrected MM=%g (corrected MM)/(original MM)=%g",MM,MM/sat); }

  arrayzero(chi,4);

  chi[0].M=factor*MM;
  chi[0].dM=factor*dMM;
  chi[0].Nm=Nm;

  sumM=sumMM=0;
  loop (i,0,3) {
    char *Mc=string(Mfmt,i+'x');

    VarM[i]=StaVar(Mc);
    if (StaError) {
      WARNING(("Variances by components: <*M%c> not found, skipped",i))
        break; }
    Nm=StaN(Mc);
    sumM+=VarM[i];
    sumMM+=Sqr(VarM[i]); }

  prt("Var M_i = %g %g %g  sum = %g",VARG(VarM),sumM);

  chi[1].M=factor*sumM;
  chi[1].dM=factor*sqrt((sumMM/3-Sqr(sumM/3))/2);
  chi[1].Nm=Nm;

  nchi=2;

#ifdef POLAR
  chi_alpha=No.A*(4*PI/V);
  put(chi_alpha)
  if (option('p')/10%10==0) {
    chi_alpha=0;
    prt("Car-Parrinello-like: chi_alpha=0 (eps_hf=1) assumed (fluctuations of\n\
*** induced dipoles are described explicitly for thermalized Drude oscilators)"); }
  chi[0].eps_hf=chi[1].eps_hf=(3+2*chi_alpha)/(3-chi_alpha);
  chi[0].hf=chi[1].hf=(2*el.epsinf+1)*(chi[0].eps_hf-1)/(chi[0].eps_hf+2*el.epsinf);

  chi[2].hf=StaMean("high-frequency chi");
  nchi+=!StaError;
  if (nchi>2) {
    chi[2].dhf=StaStdErr("high-frequency chi");
    chi[2].M=chi[0].M;
    chi[2].dM=chi[0].dM;
    chi[2].Nm=StaN("high-frequency chi"); }

  chi[3].hf=StaMean("total chi");

  nchi+=!StaError;
  if (nchi>3) {
    double x;

    if (Eext.isE) {
      loop (i,0,3) if (Eext.E[i]) {
        if (scf.Estride>1) x=StaMean(string("xchi_%c",'x'+i));
        else x=sqrt(factor)*StaMean(string("M%c",'x'+i));
        put2(i,x)
        if (StaError) {
          static int pass=1;
          if (pass) WARNING(("internal: neither xchi_%c nor M%c found",'x'+i,'x'+i))
          pass=0;
          break; }
        chi[3].hf-=Sqr(x); } }
    chi[3].dhf=StaStdErr("total chi");
    chi[3].Nm=StaN("total chi");
    chi[3].M=0;
    chi[3].dM=0; }
#endif /*# POLAR */

  prt("NB: errors are expressed as:\n\
  MEAN STDERR  = values with linear estimates of standard errors\n\
  [ LOW HIGH ] = range of given confidence limit");

  loop (ichi,0,nchi) {
    double X=chi[ichi].M+chi[ichi].hf;
    double d=sqrt(Sqr(chi[ichi].dM)+Sqr(chi[ichi].dhf));
    double eps=chi2eps(X);
    double deps=(chi2eps(X+d*0.005)-chi2eps(X-d*0.005))*100;
    double eps_hf=chi2eps(chi[ichi].hf);
    double deps_hf=chi2eps(chi[ichi].hf+chi[ichi].dhf/2)-chi2eps(chi[ichi].hf-chi[ichi].dhf/2);

    switch (ichi) {
      case 0: prt("NONPOLAR or POLAR with Clausius-Mossotti-based high-frequency epsilon:"); break;
      case 1: prt("as above, error determined from three components (Var M_i):"); break;
#ifdef POLAR
      case 2:
        prt("POLAR with separately determined high-frequency epsilon\n\
* sampling elst field E = %g V/m in directions [x,-x,y,-y,z,-z]\n\
* sampling rate (scf.stride) = %d cycles = %g ps:",scf.E,scf.Estride,scf.Estride*No.dt);
        break;
      case 3:
        prt("POLAR: high-frequency response and Var M-based response measured at once:");
        if (Eext.isE) prt("WARNING: <M>^2 caused by external elst field subtracted");
        if (scf.Estride>1) prt("WARNING: low sampling rate (scf.Estride=%d)",scf.Estride);
        break;
#endif /*# POLAR */
      default: WARNING(("bad ichi - POLAR simulation continued by nonpolar cook*?")) }

    prt("  EPS%d  chi.M = %g %g  chi.hf = %g %g  chi.Nm = %d",
        ichi,chi[ichi].M,chi[ichi].dM,chi[ichi].hf,chi[ichi].dhf,chi[ichi].Nm);

    if (ichi<3) prt("  EPS%d  high-frequency epsilon = %8.6f %8.6f   (chi = %8.6f %8.6f)",
                    ichi, eps_hf,deps_hf, chi[ichi].hf,chi[ichi].dhf);
    prt("  EPS%d  total epsilon = %.4f %.4f",
        ichi,eps,deps);
    prt("  EPS%d  68%% confidence range = [ %.4f %.4f ]",
        ichi,chi2eps(X-d),chi2eps(X+d));
    prt("  EPS%d  95%% confidence range = [ %.4f %.4f ]",
        ichi,chi2eps(X-1.96*d),chi2eps(X+1.96*d));
  }

}

#ifdef SHEAR
void prtviscosity(void) /************************************** prtviscosity */
{
  if (shear) {
    double den=StaMean("shear:den");
    double dden=StaStdErr("shear:den");

    if (!StaError) {
      double num=StaMean("shear:num");
      double dnum=StaStdErr("shear:num");
      double eta=num/den;
      double relerr=sqrt(Sqr(dnum/num)+Sqr(dden/den));
      prt("\nviscosity = <shear:num>/<shear:den> = %g prog.units, rel.err=%g", eta,relerr);
      prt("viscosity/Pa.s = %g %g",eta*(Punit*timeunit),eta*(Punit*timeunit)*relerr);
    }
  }
}
#endif /*# SHEAR */

#ifdef POLAR
/* because polarrof(mn,cfg[1]->rp) is NOT velocity */
vector *lastrpols;
#endif /*# POLAR */

#ifdef RGYR
double measureplus(double *endend,int end[2],int I) /*********** measureplus */
#else /*# RGYR */
void measureplus(void) /**************************************** measureplus */
#endif /*#!RGYR */
/*
  Calculated every cycle irrespective of Ewald summation (COULOMB<0)
  - molecular dipole moments for all species separately,
    for ions with respect to the center of charge^2;
    with el.bg without background terms
  - electric current autocorrelation function (if lag.J)
  - velocity autocorrelation function (see lag.v, lag.dim)
  #ifdef RGYR also:
    gyration tensor = SUM m_i (r_i-r_CM)(r_i-r_CM)/SUM m_i
    radius of gyration = sqrt(SUM m_i (r_i-r_CM)^2/SUM m_i)
      (=return value; in statistics squared)
    |rhead-rtail|
      (=return value; in statistics squared)
*/
{
  molecule_t *mn;
  siteinfo_t *si;
  int ns,sp,gr,lastsp=-1,moln=0,n,i,sumn;
  vector *r,*v,D,CQ,CM;
  double q,sumq,sumqq,sumpolqq,sumM,sq;
  static struct sum_s {
    double D,DD; /* dipole moment per molecule, abs. and squared */
    double endend; /* end-end distance */
    vector J; /* current (by groups) */
    vector M; /* dipole moment (by groups) */
#ifdef RGYR
    double Rgyrq; /* trace of the gyration tensor = R_gyr^2 */
    vector RGm; /* eigenvalues (diagonalized) */
#endif /*# RGYR */
#ifdef POLAR
    double Dpol,DDpol,Dboth,DDboth,Dprj,DDprj;
#endif /*# POLAR */
    int nD,nRG,nend;
  } *sum; /* [ngroups] */
  static int ngroups0;
#ifdef RGYR
  double RG[3][3];
  double RGret=0,endendret=0,Rgyrqs;
  int sumnend;
  static int warning=10;
  vector *velocity;
  static double **RGa=NULL;
#endif /*# RGYR */

#ifdef POLAR
  vector Dpol,vpol;
  double qpol;
  vector *rpol,*lastrpol;
  vector sumD3,sumD3pol;
  static int warn=1;
  if (!lastrpols && warn) {
    warn=0;
    WARNING(("lastrpols not defined\n\
*** displacement current of Drude charges is not included in current density J\n\
*** - acceptable for steady state as for NEMD conductivity from external field\n\
*** - no problem when called from plb2diff and similar because not used at all")) }
#  define sumD3both En.M
  VVO(sumD3pol,=sumD3both,=0)
#else /*# POLAR */
#  define sumD3 En.M
#endif /*#!POLAR */

  el.sumM=0;

  /* to avoid frequent alloc/free */
  if (!ngroups0) {
    allocarrayzero(sum,ngroups);
    ngroups0=ngroups; }
  else {
    if (ngroups!=ngroups0) ERROR(("ngroups changed (old=%d,new=%d)",ngroups0,ngroups))
    arrayzero(sum,ngroups); }

#ifdef RGYR
  if (!RGa) alloc2Darray(RGa,3,3);
  Min(lag.nv,No.N)
  if (lag.dim>DIM) {
    WARNING(("lag.dim=%d out of range: changed to 3",lag.dim))
    lag.dim=3; }
  if (lag.v) allocarrayzero(velocity,lag.nv);
#endif /*# RGYR */

  VO(sumD3,=0)

  loop (n,0,No.N) {
    mn=molec+n;
    ns=mn->ns;
    sp=mn->sp;
    if (sp==lastsp) moln++;
    else moln=0;
    lastsp=sp;
    si=spec[sp]->si;
    gr=spec[sp]->group;
    r=rof(mn,cfg[0]->rp);
    v=rof(mn,cfg[1]->rp);
#ifdef POLAR
    rpol=polarrof(mn,cfg[0]->rp);
    /* polarrof(mn,cfg[1]->rp) is NOT velocity */
    if (lastrpols) lastrpol=((vector*)((char*)(lastrpols)+mn->ir));
    VO(Dpol,=0)
#endif /*# POLAR */
    VVVO(D,=CQ,=CM,=0)
    sumq=sumqq=sumpolqq=sumM=0;
#ifdef RGYR
    memset(RG,0,sizeof(RG));
#endif /*# RGYR */
    loop (i,0,ns) {
      q=si[i].charge;
      /* dipole moments */
#ifdef POLAR
#  if POLAR&32
      if (si[i].qtype&32) {
        /* to be doublechecked: */
        qpol=rpol[i][0]; /* q + delta q */
        if (el.centroid) VV(CQ,+=Sqr(qpol)*r[i]) /* ions: sqr charge centroid */
        sumqq+=Sqr(qpol);
        sumq+=qpol; /* NB: sum qpol = sum q */
        VV(sum[gr].J,+=qpol*v[i])
        VV(sum[gr].M,+=qpol*r[i])
        VV(D,+=q*r[i])
        VV(Dpol,+=(qpol-q)*r[i]) }
      else
      /* optimization possible - switch nonpol, skip uncharged... */
#  endif /*# POLAR&32 */
        { /* Drude */
          qpol=si[i].chargepol;

          /* current density */
          if (lastrpols) {
            /* exact incl. the displacement current */
            VVV(vpol,=rpol[i],-lastrpol[i])
            VVV(sum[gr].J,+=(q+qpol)*v[i],+qpol*vpol) }
          else
            /* approximate excl. the displacement current (good enough for NEMD) */
            VV(sum[gr].J,+=(q+qpol)*v[i])

          /* permanent dipoles */
          q+=si[i].chargepol;
          if (el.centroid) VV(CQ,+=Sqr(q)*r[i]) /* ions: sqr charge centroid (w/o Drude) */
          sumqq+=Sqr(q);
          sumpolqq+=Sqr(qpol);
          sumq+=q;
          VV(D,+=q*r[i])
          VV(sum[gr].M,+=q*r[i])

          /* induced dipoles */
          VV(Dpol,+=qpol*rpol[i])
          VV(sum[gr].M,+=qpol*rpol[i])
        }
#else /*# POLAR */
        /* nonpolar */
        if (el.centroid) VV(CQ,+=Sqr(q)*r[i]) /* ions: sqr charge centroid */
        sumqq+=Sqr(q);
        VV(sum[gr].J,+=q*v[i])
        VV(sum[gr].M,+=q*r[i])
        sumq+=q;
        VV(D,+=q*r[i])
#endif /*#!POLAR */

#ifdef RGYR
      /* summing up CM+RG things */
      q=si[i].mass;
      sumM+=q;
      VV(CM,+=q*r[i])
      if (lag.v) VV(velocity[n],+=q*v[i])
      {
        int ii,jj;
        /* NB: symmetry used */
        loop (ii,0,3) loopto (jj,0,ii) RG[ii][jj]+=q*r[i][ii]*r[i][jj];
      }
#endif /*# RGYR */
    } /* i */

    if (sumqq+sumpolqq) {
      sum[gr].nD++;
      /* dipole moment for ions to q^2 centroid (w/o Drude charge) */
      if (el.centroid && sumqq!=0) VV(D,-=sumq/sumqq*CQ)

      /* permanent dipole (only this one for nonpolarizable models) */
      q=SQR(D);
      sum[gr].D+=sq=sqrt(q);
      if (sq==0) sq=1; /* to avoid nan */
      sum[gr].DD+=q;
      VV(sumD3,+=D)
#ifdef POLAR
      /* induced dipole */
      qpol=SQR(Dpol);
      sum[gr].Dpol+=sqrt(qpol);
      sum[gr].DDpol+=qpol;
      VV(sumD3pol,+=Dpol)

      /* projection of the induced dipole to the direction of the permanent one */
      qpol=SCAL(D,Dpol)/sq;
      sum[gr].Dprj+=qpol; /* bug fixed */
      sum[gr].DDprj+=Sqr(qpol);

      /* both permanent + induced */
      VV(D,+=Dpol)
      q=SQR(D);
      sum[gr].Dboth+=sqrt(q);
      sum[gr].DDboth+=q;
      VV(sumD3both,+=D)
#endif /*# POLAR */
    } /* sumqq */

#ifdef RGYR
    if (lag.v && n<lag.nv) {
      StaSet(0,lag.v,2,0);
      loop (i,0,lag.dim)
        StaAdd(string("v%c%d.%d",i+'x',sp,moln),velocity[n][i]/sumM); }

    if (ns>1) {
      double RGval,endendval,x;
      static int End[2],oldsp=-1;
      int ii,jj;

      /* filling the symmetric parts */
      RG[1][2]=RG[2][1];
      RG[0][1]=RG[1][0];
      RG[0][2]=RG[2][0];

      VO(CM,/=sumM)
      RGval=0;
      loop (ii,0,3) {
        loop (jj,0,3)
          RGa[ii][jj]=RG[ii][jj]/sumM - CM[ii]*CM[jj];
        RGval+=RGa[ii][ii]; }
      sum[gr].Rgyrq += RGval;

      Jacobi(-3,RGa,NULL,No.eps);

      /* order Rgxx<Rgyy<Rgzz */
      if (RGa[0][0]>RGa[1][1]) x=RGa[0][0],RGa[0][0]=RGa[1][1],RGa[1][1]=x;
      if (RGa[1][1]>RGa[2][2]) x=RGa[1][1],RGa[1][1]=RGa[2][2],RGa[2][2]=x;
      if (RGa[0][0]>RGa[1][1]) x=RGa[0][0],RGa[0][0]=RGa[1][1],RGa[1][1]=x;

      loop (ii,0,3) sum[gr].RGm[ii]=RGa[ii][ii];

      sum[gr].nRG++;

      if (sp!=oldsp) {
        if (sp<oldsp) warning=0;
        loop (i,0,2) {
          if (warning && (end[i]<-ns || end[i]>=ns)) {
            prt("WARNING: molecule %d, ns=%d: bad rg.end[%d]=%d",
                n,ns,i,end[i]);
            if (--warning==0) prt("... more warnings suppressed"); }

          if (end[i]<0) End[i]=ns+end[i]; else End[i]=end[i];
          if (End[i]<0) End[i]=0;
          if (End[i]>=ns) End[i]=ns-1; }

        if (End[0]==End[1]) End[0]=0,End[1]=ns-1;
        oldsp=sp; }

      VVV(CM,=r[End[0]],-r[End[1]])
      sum[gr].nend++;
      sum[gr].endend += (endendval=SQR(CM));

      if (I==n) RGret=RGval,endendret=endendval; }
    else
      sum[gr].endend=RGret=endendret=0;
#endif /*# RGYR */
  } /* n */

#ifdef RGYR
  warning=999;
#endif /*# RGYR */

  /* by groups: J=el. current density (not current) */
  VO(En.J,=0)
  q=(chargeunit/Sqr(lengthunit)/timeunit)/(PROD(box.L)*h);
  loop (n,0,ngroups) {
    VO(sum[n].J,*=q)
    if (n<3) VV(En.Jg[n],=sum[n].J);
    if (lag.J>0) {
      StaSet(0,lag.J,2,(lag.J>1)*lag.n);
      StaAdd(string("Jx[%d] [A/m2]",n),sum[n].J[0]);
      StaAdd(string("Jy[%d] [A/m2]",n),sum[n].J[1]);
      StaAdd(string("Jz[%d] [A/m2]",n),sum[n].J[2]); }
    VV(En.J,+=sum[n].J) }

  if (lag.J) {
    StaSet(0,abs(lag.J),2,(abs(lag.J)>1)*lag.n);
    StaAdd("Jx [A/m2]",En.J[0]);
    StaAdd("Jy [A/m2]",En.J[1]);
    StaAdd("Jz [A/m2]",En.J[2]);
    StaAdd("J^2 [A2/m4]",SQR(En.J)); }

  /* by groups: M=dipole moment, might be with jumps for ions */
  q=1; /* p.u. kept - change into polarization? */
  loop (n,0,ngroups) {
    VO(sum[n].M,*=q)
    if (n<3) VV(En.Mg[n],=sum[n].M);
    if (lag.M>0) {
      StaSet(0,lag.M,2,(lag.M>1)*lag.n);
      StaAdd(string("Mx[%d]",n),sum[n].M[0]);
      StaAdd(string("My[%d]",n),sum[n].M[1]);
      StaAdd(string("Mz[%d]",n),sum[n].M[2]); } }

  if (lag.M) {
    StaSet(0,abs(lag.M),2,(abs(lag.M)>1)*lag.n);
    StaAdd("Mx",En.M[0]);
    StaAdd("My",En.M[1]);
    StaAdd("Mz",En.M[2]);
    StaAdd("M^2",SQR(En.M)); }

  sumn=0;
  StaSet(0,2,2,lag.n);
#ifdef POLAR
  loop (n,0,ngroups) if (sum[n].nD) {
    sumn++;
    StaAdd(string("mu_perm[%d]",n),sum[n].D/sum[n].nD);
    /* sumM does not include saturation: el.sumM+=sum[n].D; */
    StaAdd(string("mu_perm[%d]^2",n),sum[n].DD/sum[n].nD);
    StaAdd(string("mu_ind[%d]",n),sum[n].Dpol/sum[n].nD);
    StaAdd(string("mu_ind[%d]^2",n),sum[n].DDpol/sum[n].nD);
    StaAdd(string("mu_tot[%d]",n),sum[n].Dboth/sum[n].nD);
    /* sumM includes saturation: bug fixed in 2.7h */
    el.sumM+=sum[n].Dboth;
    StaAdd(string("mu_tot[%d]^2",n),sum[n].DDboth/sum[n].nD);
    StaAdd(string("mu_prj[%d]",n),sum[n].Dprj/sum[n].nD);
    StaAdd(string("mu_prj[%d]^2",n),sum[n].DDprj/sum[n].nD); }

  if (sumn) {
    StaAdd("(sum mu_perm)^2",SQR(sumD3));
    StaAdd("(sum mu_ind)^2",SQR(sumD3pol));
    StaAdd("(sum mu_tot)^2",SQR(sumD3both)); }
#else /*# POLAR */
  loop (n,0,ngroups) if (sum[n].nD) {
    sumn++;

    StaAdd(string("mu[%d]",n),sum[n].D/sum[n].nD); el.sumM+=sum[n].D;
    StaAdd(string("mu[%d]^2",n),sum[n].DD/sum[n].nD); }

  if (sumn) StaAdd("(sum mu)^2",SQR(sumD3));

  /* magnetic dipole moment */
  if (Eext.ism) {
    vector magm={0,0,0},dr;
    double rr,qm,maxmag=0;

    loop (n,0,No.N) {
      mn=molec+n;
      r=rof(mn,cfg[0]->rp);
      if (mn->sp==el.m.sp || (el.m.sp<0 && el.m.sp+mn->sp>=0)) {
        VVV(dr,=r[el.m.plus],-r[el.m.minus])
        rr=sqrt(SQR(dr));
        qm=Eext.m/rr;
        maxmag+=Eext.m;
        VV(magm,+=qm*dr) } }
    StaAdd("maxmagnetization [A/m]",maxmag/box.V*(chargeunit/lengthunit/timeunit));
    loop (n,0,DIM)
      StaAdd(string("magnetization[%d] [A/m]",n),magm[n]/box.V*(chargeunit/lengthunit/timeunit));
  }
#endif /*#!POLAR */

#ifdef RGYR
  sumn=sumnend=0;
  VO(En.RGm,=0)
  Rgyrqs=*endend=0;
  loop (n,0,ngroups) {
    if (sum[n].nRG) {
      Rgyrqs+=sum[n].Rgyrq;
      VV(En.RGm,+=sum[n].RGm)
      sumn+=sum[n].nRG;
      StaAdd(string("Rgxx[%d]",n),sum[n].RGm[0]/sum[n].nRG);
      StaAdd(string("Rgyy[%d]",n),sum[n].RGm[1]/sum[n].nRG);
      StaAdd(string("Rgzz[%d]",n),sum[n].RGm[2]/sum[n].nRG);
      StaAdd(string("Rgyr[%d]^2",n),sum[n].Rgyrq/sum[n].nRG); }
    if (sum[n].nend) {
      *endend+=sum[n].endend;
      sumnend+=sum[n].nend;
      StaAdd(string("endend[%d]^2",n),sum[n].endend/sum[n].nend); } }
  VO(En.RGm,/=sumn)

  if (lag.v) free(velocity);

  if (RGret>0) {
    *endend=sqrt(endendret);
    return sqrt(RGret); }
  else {
    if (sumnend==0) *endend=0; else *endend=sqrt(*endend/sumnend);
    if (sumn==0) return 0; else return sqrt(Rgyrqs/sumn); }
#endif /*# RGYR */
}

void ssdistance(void) /*************************************** distancecheck */
{
  if (ssd0) {

    struct ssd_s *ssd;
    vector *r=cfg[0]->rp,dr;

    StaSet(0,2,2,lag.n);

    for (ssd=ssd0; ssd; ssd=ssd->next) if (ssd->indx[0]>=0) {
      VVV(dr,=r[ssd->indx[0]],-r[ssd->indx[1]])
#ifndef FREEBC
      {
        int k;
        loop (k,0,3) {
          while (dr[k]>box.Lh[k]) dr[k]-=box.L[k];
          while (dr[k]<-box.Lh[k]) dr[k]+=box.L[k]; }
      }
#endif /*# FREEBC */
      ssd->u.dist=sqrt(SQR(dr));
      /* postponed: StaAdd(ssd->name,ssd->u.dist); */ } }
}

#ifdef HARMONICS
#  if HARMONICS==3
#    include "harmwat.c"
#  elif HARMONICS==2
#    include "harmlin.c"
#  elif HARMONICS==1
#    include "harmmeoh.c"
#  else /*#!HARMONICS==3!HARMONICS==2!HARMONICS==1 */
#    error "wrong value of HARMONICS: use 3=water, 2=linear molecule"
#  endif /*#!HARMONICS==3!HARMONICS==2!HARMONICS==1 */
#endif /*# HARMONICS */

/****************************************************************************/

#ifdef WIDOM
/*
  WIDOM insertion and scaled particle (growing molecule or identity change).
  Note that the LINKCELL method calculates the insertion a and growth
  energies using the configuration BEFORE the current integration step is
  performed while the standard version AFTER it.
*/

#  ifdef POLAR
#    warning "WIDOM+POLAR is not a good idea (unless uncharged molecule is inserted)"
#  endif /*# POLAR */
/* code doubled -- unify w. simils !!! */
/* USE HEADERS!!! */
double rndcos();
double rnd();

extern vector **initsites;

#  ifdef SLAB
double mol_wall(vector f[],molecule_t *m,vector *rp);

double Widom(int sp,int ncyc,double z0,double z1,double dz,int mode)
/*                                                                   * Widom *
   mode&1: incl. symmetrization
   mode&2: total chem.pot.
   mode&4: residual (without wall potential)
   mode&8: wall potential (angle-averaged)
*/
#  else /*# SLAB */
double Widom(int sp,int ncyc,double corr) /*************************** Widom */
#  endif /*#!SLAB */
{
  En_t En0=En;
  rdf_t ***rdf0=rdf;
  int icyc,ns=spec[sp]->ns,i;
  static vector *rp,*rpf,*c;
  static molecule_t mol;
  double ps[3];
  double ret=0;
  vector dr;

#  ifdef SLAB
  double z,Ewall=0;
  double (*w)[3];
  int nz=(z1-z0)/dz+0.9999,iz,k;

  if (nz<=0 || ncyc<=0) return 0;
  if (mode&1 && fabs((z1-z0)/dz-nz)>1e-4)
    ERROR(("Widom: z1-z0 must be an integer multiple of dz if mode&1"))

  if (mode&1) allocarray(w,nz+1);
#  else /*# SLAB */
  /* cutoff correction - multiplication factor
     it is V-dependent, thus NPT is allowed though not recommended */
  double qcorr=exp(-corr/(box.V*T));

  if (ncyc<=0) return 0;
#  endif /*#!SLAB */

  rdf=NULL; widomrdf(0);

  if (!rp) {
    allocarray(rp,No.s+ns); /* cfg+probe */
    allocarrayzero(rpf,No.s+ns); /* dummy forces */ }
  mol.sp=sp;
  mol.ns=ns;
  mol.nc=spec[sp]->nc;
  mol.ir=No.s*sizeof(vector);
  c=rp+No.s;
  mol.ig=0; /* not needed */

  StaSet(0,2,2,lag.n);

#  ifdef SLAB
  loopto (iz,0,nz) {
    z=z0+iz*dz;
#  endif /*# SLAB */

    ps[0]=ps[1]=ps[2]=0;

#  ifndef LINKCELL
    depend_r(cfg[0]); /* why this? */
    memcpy(rp,cfg[0]->rp,No.s*sizeof(vector));
#  endif /*# LINKCELL */

#  ifdef SLAB
    loop (icyc,0,ncyc) {
      dr[0]=rnd()*L[0];
      dr[1]=rnd()*L[1];
      dr[2]=z;
#  else /*# SLAB */
    loop (icyc,0,ncyc) {
      VV(dr,=rnd()*box.L)
#  endif /*#!SLAB */

      loop (i,0,ns) VV(c[i],=initsites[sp][i])
      rndorientation(c,ns);

      loop (i,0,ns) VV(c[i],+=dr)

#  ifdef LINKCELL
      if (ns>1) {
        /* the molecule position must be exactly normalized, i.e.,
           no coordinate may be negative; cf. normalize(int m) in norm.c
           (the direct (all-pairs) method is not sensitive
           to small negative cordinates)
           (it is assumed that 1 atom is always centered!!!) */
        vector rmin;
        int k;

        VV(rmin,=c[0])
        loop (i,1,ns) loop (k,0,DIM)
          if (c[i][k]<rmin[k]) rmin[k]=c[i][k];
        loop (k,0,DIM)
          if (rmin[k]<0) loop (i,0,ns) c[i][k]+=box.L[k]; }
#  endif /*# LINKCELL */

#  ifdef SLAB
      if (wall.is) {
        En.el=0;
        Ewall=mol_wall(rof(&mol,rpf),&mol,rp);
        Ewall=-(Ewall+En.el)/T; }
#  endif /*# SLAB */

      /* insertion energy returned in En.pot */
      Einsert(&mol,rpf,rp);
      En.pot/=-T; /* WARNING: Widom requires thermostat! */

#  ifdef SLAB
      if (mode&2) if (En.pot>-700) ps[0]+=exp(En.pot);
      if (mode&4) if (Ewall+En.pot>-700) ps[1]+=exp(En.pot+Ewall);
      if (mode&8) if (Ewall>-700) ps[2]+=exp(Ewall); }

    if (mode&4) ret+=ps[1];
    else if (mode&2) ret+=ps[0];
    else ret+=ps[2];

    loop (k,0,3) {
      int m=2<<k;
      ps[k]/=ncyc;
      if (mode&m) StaAdd(string("Widom%d z=%g",m,z),ps[k]);
      if (mode&1) w[iz][k]=ps[k]; } } /* iz-loop */
    ret/=ncyc*nz;
#  else /*# SLAB */
    if (En.pot>-700) ps[0]+=exp(En.pot); }
  StaAdd("Widom raw",ps[0]/ncyc);
  StaAdd("Widom",ps[0]/ncyc*qcorr);
  ret=ps[0]/ncyc*qcorr;
#  endif /*#!SLAB */

#  ifdef SLAB
  if (mode&1) {
    loop (iz,0,nz/2) {
      z=z0+iz*dz;
      loop (k,0,3) {
        int m=2<<k;
        if (mode&m)
          StaAdd(string("Widom%d z=%g",m+1,z),(w[iz][k]+w[nz-iz][k])/2); } }
    free(w); }
#  endif /*# SLAB */

  En=En0;
  rdf=rdf0;
  widomrdf(1);

  return ret;
}


double XWidom(int sp,int spvirt) /********************************** XWidom */
{
  En_t En0=En;
  rdf_t ***rdf0=rdf;
  int n,nscale=0;
  static vector *rp,*rpf;
  double ret=0;
  double e,psum=0;

  rdf=NULL; widomrdf(0);

  if (!rp) {
    allocarray(rp,No.s); /* cfg */
    allocarrayzero(rpf,No.s); /* dummy */ }

  StaSet(0,2,2,lag.n);

#  ifndef LINKCELL
  depend_r(cfg[0]); /* why this? */
  memcpy(rp,cfg[0]->rp,No.s*sizeof(vector));
#  endif /*# LINKCELL */

  loop (n,0,No.N) if (molec[n].sp==sp) {
    nscale++;
    Eidchange(spvirt,&molec[n],rpf,rp);
    En.pot/=-T; /* WARNING: Widom requires thermostat! */
    e=0; if (En.pot>-700) e=exp(En.pot);
#  ifdef SLAB
    StaAdd(string("XWidom%d",n),e);
#  endif /*# SLAB */
    psum+=e; }

  ret=psum/nscale;
#  ifndef SLAB
  StaAdd("XWidom",ret);
#  endif /*# SLAB */

  En=En0;
  rdf=rdf0;
  widomrdf(1);

  return ret;
}
#endif /*# WIDOM */

#ifdef SLAB
#  include "slabmeas.c"
#endif /*# SLAB */

// #ifdef SLAB
// #  if SLAB & 8
// #    include "slabmeas-mapaver.c"
// #  else /*# SLAB & 8 */
// #    include "slabmeas.c"
// #  endif /*#!SLAB & 8 */
// #endif /*# SLAB */
