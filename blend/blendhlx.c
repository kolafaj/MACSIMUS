/* alpha-helix + backbone using phi,psi,omega */

void alphahelix(species_t *spec) /******************************* alphahelix */
{
 int i,ns=spec->ns,nbb=0,nsc=0,j;
 site_t *si,*sj;
 vector r;

 double alpha=0,z=0,rad,al;
 const double dalpha=0.582; /* angle/bead = 3.6 residues/turn */
 const double dz=0.5;       /* delta z/bead in A =1.5/residue */
 const double radius=2.3;   /* radius in A (of C-alpha) */

 makersdinfo(spec);

 spec->opt_k=1;
 
 VO(r,=0)

 site=spec->site;

 loop (i,0,ns)
   if ( (si=spec->site+i)->backbone && si->backbone<8) {
     si->keep=INJAIL;
     rad=radius;
     r[2]=z;
     al=alpha;
     /* N=1 Calpha=2 C=3 */
     if (si->backbone==1 || si->backbone==3) {
       j=i;
       rad*=0.72;
       if (si->backbone==1) r[2]-=0.5; else r[2]+=0.5; }
     if (si->backbone==2) al-=.2;
     if (si->backbone==3) al-=.3;
     r[0]=cos(alpha)*rad;
     r[1]=sin(alpha)*rad;
     z+=dz; alpha+=dalpha;
     copy(si->r,r,sizeof(vector));
     if (si->backbone==2) {
/*.....       loop (j,0,ns) { sj=spec->site+j; prt("%d %d %d",sj->rsd,sj->backbone,sj->nnbr); }*/
       loop (j,0,ns) if ( (sj=spec->site+j)->rsd==si->rsd && !sj->backbone
                          && sj->nnbr>1) {
         /* very stupid: just one chain, bad for rings */
         sj->keep=FREE;
         sj->r[0]=r[0]*(1+sj->nest+(j-i)*.01)*.5;
         sj->r[1]=r[1]*(1+sj->nest-(j-i)*.01)*.5;
         sj->r[2]=r[2]-1; 
         VO(sj->r,+=0.05*rndcos()) }
       nsc++; }
     nbb++; }
   else
     si->keep=WANTED;


#if 0
 loop (i,0,ns) 
   prt("%3d %4s bb=%2d rsd=%2d keep=%d",
       i,atom[spec->site[i].type],spec->site[i].backbone,spec->site[i].rsd,spec->site[i].keep);
#endif

 prt("! alpha-helix initial configuration, %d in backbone, %d in side-chains",
     nbb,nsc);
} /* alphahelix */

void alphahelixold(species_t *spec) /************************* alphahelixold */
/* ... not very sophisticated.  Assumes CHARMM21 atom names. */
/* NOTE: to be updated using function backbone() */
{
int t,i,j=-1,ns=spec->ns,nbb=0,nsc=0;
site_t *si;
vector r;

double alpha=0,z=0,rad;
const double dalpha=0.582; /* angle/bead = 3.6 residues/turn */
const double dz=0.5;       /* delta z/bead in A =1.5/residue */
const double radius=2.5;   /* radius in A (of C-alpha) */

spec->opt_k=1;

VO(r,=0)

loop (i,0,ns)
  if ((si=spec->site+i)->nnbr<2)
    si->keep=WANTED;
  else {
    if ( (t=bbtype(i)) ) if (isnbr(j,si)) {
      j=i;
      si->keep=INJAIL;
      rad=radius;
      r[2]=z;
      if (t!=2) {
        rad*=0.72;
        if (t==1) r[2]-=0.5; else r[2]+=0.5; }
      r[0]=cos(alpha)*rad;
      r[1]=sin(alpha)*rad;
/* prt("%d %2d %4s : %.3f",t,si->type,typenames[t],r[2]); */
      z+=dz; alpha+=dalpha;
      copy(si->r,r,sizeof(vector));
      nbb++;
      goto done; }

    si->keep=FREE;
    /* very stupid: just one chain, bad for rings */
    si->r[0]=r[0]*(1+i-j);
    si->r[1]=r[1]*(1+i-j);
    si->r[2]=r[2]-1;
    VO(si->r,+=0.05*rndcos())
    nsc++;
    done:; }

prt("!OLD: alpha-helix initial configuration, %d in backbone, %d in side-chains",
    nbb,nsc);
} /* alphahelixold */

static vector O[3]; /* orientation */
static void rotate(int x,double angle)
{
  vector OO[3];
  int y=(x+1)%3,z=(y+1)%3;
  double ca,sa;

  angle*=PI/180;
  ca=cos(angle),sa=sin(angle);

  copy(OO,O,sizeof(OO));
  O[y][x]=OO[y][x]*ca+OO[z][x]*sa;
  O[y][y]=OO[y][y]*ca+OO[z][y]*sa;
  O[y][z]=OO[y][z]*ca+OO[z][z]*sa;
  O[z][x]=OO[z][x]*ca-OO[y][x]*sa;
  O[z][y]=OO[z][y]*ca-OO[y][y]*sa;
  O[z][z]=OO[z][z]*ca-OO[y][z]*sa;
  /*
  prt("%9.3f %9.3f %9.3f",O[0][0],O[0][1],O[0][2]);
  prt("%9.3f %9.3f %9.3f",O[1][0],O[1][1],O[1][2]);
  prt("%9.3f %9.3f %9.3f",O[2][0],O[2][1],O[2][2]);
  _n
  */
}

void makechain(species_t *spec) /********************************* makechain */
{
  int i,ns=spec->ns,nbb=0,nsc=0,x,y;
  site_t *si;
  vector R; /* translation */
  static ireal bbond[3]={1.46,1.515,1.335}; /* N-Ca, Ca-C, C-N */
  static ireal bangle[3]={180-112.8, /* N-Ca-C */
			  180-113.9, /* Ca-C-N */
			  180-117.6};/* C-N-Ca */
  makersdinfo(spec);

  spec->opt_k=1;
 
  VO(R,=0)
  loop (x,0,3) loop (y,0,3) O[x][y]=x==y;

  site=spec->site;

  loop (i,0,ns) spec->site[i].keep=WANTED;

  loop (i,0,ns)
    if ( (si=spec->site+i)->backbone>0 && si->backbone<4) {
     /* N=1 Calpha=2 C=3: it is asssumed that sites come in this order! */
      si->keep=INJAIL;
      VV(si->r,=R)
      rotate(1,bangle[si->backbone-1]);
      VV(R,+=bbond[si->backbone-1]*O[2])
      rotate(2,Xopt.bdihedral[si->backbone-1]);
      
      nbb++; }

  prt("! chain initial configuration, %d in backbone, %d in side-chains",
      nbb,nsc);
} 

