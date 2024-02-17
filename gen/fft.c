/*
  Generalized fast Fourier transform
  History :
    2016 small improvement if nh=product of large primes
         PARALLEL (via #pragma parallel etc.) removed
    2008 C-module, 3D complex FFT added
    cca. 1995 fourier.c (now deprecated)
    1990 FOURIER.PAS

Declarations required:
^^^^^^^^^^^^^^^^^^^^^^
  Type  fftreal == <real | single | double | extended>;
  Const Max3Dh == <max grid for (3D>; ($IfDef ER3DL || ER3D only)

Procedures available:
^^^^^^^^^^^^^^^^^^^^^
  Fourier(var a:CVectorPtr; n:integer);
  FourierER(a:RVectorPtr; nh,by:integer);        ($IfDef ER )
  FourierER4(a,b,c,d:RVectorPtr; nh,by:integer); ($IfDef ER )
  FourierER3D(var a:VectorER3D; nh:integer);     ($IfDef ER && ER3D)
  FourierER3DL(var a:Large3D; nh:integer);       ($IfDef ER && ER3DL)

============================ FOURIER.PAS ==============================*/

#include "ground.h"
#include "fft.h"

QPtr QHead=NULL;

QPtr InitFour(int n) /*============================================ InitFour */
/* returns pointer to the table of Q[j]==exp(i*j/n), j==0..n-1,
  creates new table if (necessary */
{
  fftreal x;
  int i;
  QPtr Qn;

  Qn=QHead; /* searching for (Q-table */
  while (Qn!=NULL)
    if (Qn->NofQ==n) return Qn;
    else Qn=Qn->next;

  /* creating new Q-table */
  alloc(Qn,sizeof(Qtype)+(abs(n)-1)*sizeof(fftcomplex));
  Qn->next=QHead; QHead=Qn;
  Qn->NofQ=n;
  loop (i,0,abs(n)) {
    x=2*PI/n; Qn->Q[i].re=cos(x*i); Qn->Q[i].im=sin(x*i); }

  return Qn;
} /* InitFour */


void Fourier(fftcomplex *a, int n) /*================================== Fourier */
/*
fftcomplex array a[0..m-1], m==abs(n), is replaced by the Fourier transform:

           m-1
    a[k] = SUM a[j]*exp[i*j*k/n], k==0..m-1
           j=1

NOTES:
* n<0 denotes back transform.
* Efficiency decreases if  n  contains large prime numbers,
  recommended values are  n == 2**k1 * 3**k2, where k2 is not too big. */
{
  int iaux,iw,ia,iq,i,j,k,l,m,d,blk,sqblk,NewBlk,NdD,jBlk,jNewBlk,kNdD;
  unsigned int bLen;
  fftcomplex aux,w,ww,w1,w2,ww1,ww2,r3;
  fftcomplex *aparm,*b,*swap,*Q;
  QPtr Qn;

  r3.re=-0.5; r3.im=sqrt(0.75);
  if (n<0) r3.im=-r3.im;
  Qn=InitFour(n);
  n=abs(n);
  blk=n;
  bLen=n*sizeof(fftcomplex); alloc(b,bLen);

/* orig:
REPEAT
  a2b(a^,b^);
  if (blk>1) a2b(b^,a^) else move(b^,a^,bLen);
  UNTIL blk<==1;
*/

  aparm=a;

  do { /* one step */

    sqblk=(int)sqrt(blk+0.5);
    loopto (d,2,sqblk) if (blk % d == 0) break;

    if (blk % d != 0) d=blk; /* blk = prime */

    m=n/blk; NewBlk=blk/d; NdD=n/d;

    jBlk=0; jNewBlk=0;
    loop (j,0,m) {
      loop (i,0,NewBlk)
        switch (d) {

          case 2: /* optimized code */
            ia=jBlk+i;
            aux=a[ia+NewBlk];
            Q = Qn->Q + jNewBlk;
            w.re = Q->re*aux.re - Q->im*aux.im;
            w.im = Q->im*aux.re + Q->re*aux.im;
            iw=jNewBlk+i;
            aux=a[ia];
            b[iw].re=aux.re+w.re; b[iw].im=aux.im+w.im;
            b[iw+NdD].re=aux.re-w.re; b[iw+NdD].im=aux.im-w.im;
            break; /*2*/

         case 3: /* optimized code */
            ia=jBlk+i; iq=jNewBlk; iaux=ia+NewBlk;
            aux=a[iaux];
            Q = Qn->Q + iq;
            w.re = Q->re*aux.re - Q->im*aux.im;
            w.im = Q->im*aux.re + Q->re*aux.im;
            aux=a[iaux+NewBlk];
            Q = Qn->Q + iq*2 % n;
            ww.re = Q->re*aux.re - Q->im*aux.im;
            ww.im = Q->im*aux.re + Q->re*aux.im;
            iw=jNewBlk+i;
            aux=a[ia];
            b[iw].re=aux.re+w.re+ww.re; b[iw].im=aux.im+w.im+ww.im;
            w1.re=w.re*r3.re; w2.re=w.im*r3.im;
            w2.im=w.re*r3.im; w1.im=w.im*r3.re;
            ww1.re=ww.re*r3.re; ww2.re=ww.im*r3.im;
            ww2.im=ww.re*r3.im; ww1.im=ww.im*r3.re;
            iaux=iw+NdD;
            b[iaux].re=aux.re+w1.re-w2.re+ww1.re+ww2.re;
            b[iaux].im=aux.im+w1.im+w2.im+ww1.im-ww2.im;
            iaux += NdD;
            b[iaux].re=aux.re+w1.re+w2.re+ww1.re-ww2.re;
            b[iaux].im=aux.im+w1.im-w2.im+ww1.im+ww2.im;
            break; /*3*/

          default: /* d>3: not optimized */
            kNdD=0; iw=jBlk+i;
            loop (k,0,d) {
              ia=iw;
              w=a[ia];
              iaux=jNewBlk+kNdD; iq=iaux;
              loop (l,1,d) {
                Q = Qn->Q + iq;
                ia += NewBlk; iq=(iq+iaux) % n;
                aux=a[ia];
                w.re += aux.re*Q->re - aux.im*Q->im;
                w.im += aux.re*Q->im + aux.im*Q->re; }
              b[jNewBlk+i+kNdD]=w;
              kNdD += NdD; } /*k*/

          } /* switch(d), loop(i) */
        jBlk += blk; jNewBlk += NewBlk;
      } /*j*/

    blk=NewBlk;

    swap=a; a=b; b=swap;
    } while (blk>1);

  /* to preserve the stack-like structure of the heap */
  if (a!=aparm) {
    swap=a; a=b; b=swap;
    copy(a,b,bLen); }

  free(b);
} /* Fourier */


void FourierER(fftreal *a, int nh, int by) /* =================== FourierER */
/*
  Fourier transform of 1 Even Real function a^ [0..nh] - not optimized !
  Period of the function is n==2*nh.
  'a^' may be a non-contiguous vector, increment is 'by' in units of
  fftreal size (by==1 corresponds to a usual contiguous vector).
*/
{
  fftcomplex *f;
  unsigned int fLen;
  int i,n;

  InitFour(nh*2);
  fLen=sizeof(fftcomplex)*nh*2; alloc(f,fLen);
  nh=abs(nh); n=2*nh;

  loopto (i,0,nh) { f[i].re=a[i*by]; f[i].im=0; }
  loopto (i,1,nh) f[n-i]=f[i];
  Fourier(f,n);
  loopto (i,0,nh) a[i*by]=f[i].re;

  free(f);
} /* FourierER */

void FourierER4(fftreal *a, fftreal *b, fftreal *c, fftreal *d, int nh, int by)
                                                  /* ============ FourierER4 */
/*
  Fourier transform of 4 Even Real functions a^,b^,c^,d^ [0..nh].
  Optimized code, see comments of FourierER.
*/
{
  int iby,i,j,n;
  unsigned int fLen;
  fftreal c0,c1,d0,d1;
  fftcomplex fi,fni,*f,*Q;
  QPtr Qn;

  nh=abs(nh); n=nh*2;
  Qn=InitFour(n);
  fLen=sizeof(fftcomplex)*n; alloc(f,fLen);

  loop (i,0,n) {
    Q = Qn->Q + i;
    j=n-i; if (j>i) j=i; j=j*by;
    f[i].re = a[j] + Q->re*c[j] - Q->im*d[j];
    f[i].im = b[j] + Q->im*c[j] + Q->re*d[j]; }

  Fourier(f,n);

  c0=(c[0]+c[nh*by])/2; c1=(c[0]-c[nh*by])/2;
  d0=(d[0]+d[nh*by])/2; d1=(d[0]-d[nh*by])/2;

  loop (i,1,nh) {
    iby=i*by;
    Q = Qn->Q + i;
    c0=c0+c[iby]; c1=c1+c[iby]*Q->re;
    d0=d0+d[iby]; d1=d1+d[iby]*Q->re; }

  c[0]=2*c0; c[by]=2*c1; d[0]=2*d0; d[by]=2*d1;

  loop (i,1,nh) {
    iby=i*by;
    fni=f[n-i]; fi=f[i];
    c[iby+by]=fi.re-fni.re+c[iby-by];
    a[iby]=(fi.re+fni.re-c[iby-by]-c[iby+by])/2;
    d[iby+by]=fi.im-fni.im+d[iby-by];
    b[iby]=(fi.im+fni.im-d[iby-by]-d[iby+by])/2; }

  a[0]=f[0].re-c[by]; b[0]=f[0].im-d[by];
  a[nh*by]=f[nh].re-c[nh*by-by]; b[nh*by]=f[nh].im-d[nh*by-by];

  free(f);
} /* FourierER4 */

/*** for FourierER3D (optimized for x<-->y<-->z symmetric functions) ***/
void Ff4(fftreal **p, int nh, int span)
{
  FourierER4(p[0],p[1],p[2],p[3],nh,span);
}

void Fflush(fftreal **p, int nh, int span, int i)
{
  switch (i) {
    case 0: return;
    case 1: FourierER(p[0],nh,span); break;
    case 2: {
      fftreal *aux;
      alloczero(aux,(span*(nh+1))*2*sizeof(fftreal));
      FourierER4(p[0],p[1],aux,aux+span*(nh+1),nh,span);
      free(aux); } break;
    case 3: {
      fftreal *aux;
      alloczero(aux,span*(nh+1)*sizeof(fftreal));
      FourierER4(p[0],p[1],p[2],aux,nh,span);
      free(aux); } }
}

void FourierER3D(fftreal *a,int nh) /* ================== FourierER3D/serial */
/*
  3D fourier transform of Even Real 3D function a[0..nh,0..nh,0..nh].
  Period is 2*nh in each coordinate.
  For a[x,y,z]=a[y,x,z]=a[y,z,x]...

  NOTE map function: a[x,y,z] = x*span^2 + y*span + z, where span=nh+1
*/
{
  int span=nh+1;
  int x,y,z,i;
  fftreal *p[4];

  InitFour(2*nh); /* to allocate Q before Aux */

#define A(X,Y,Z) a[X*(span*span)+Y*span+Z]

  i=0;
  loopto (y,0,nh) 
    loopto (z,0,y) {
      p[i]=a+y*span+z; i++;
      if (i==4) { Ff4(p,nh,Sqr(span)); i=0; } }
  Fflush(p,nh,Sqr(span),i);
  loopto (y,0,nh) 
    loopto (z,0,y)
      loopto (x,0,nh) A(x,z,y)=A(x,y,z);

  i=0;
  loopto (x,0,nh)
    loopto (z,0,nh) {
      p[i]=a+x*Sqr(span)+z; i++;
      if (i==4) { Ff4(p,nh,span); i=0; } }
  Fflush(p,nh,span,i);

  i=0;
  loopto (x,0,nh)
    loopto (y,0,x) {
      p[i]=a+x*Sqr(span)+y*span; i++;
      if (i==4) { Ff4(p,nh,1); i=0; } }
  Fflush(p,nh,1,i);

  loopto (x,0,nh)
    loop (y,0,x)
      loopto (z,0,nh) A(y,x,z)=A(x,y,z);

#undef A
} /*FourierER3D*/


/* NEW 2008: */
void Fourier3D(fftcomplex ***A,int nx,int ny,int nz)
/*
  3D general (complex) Fourier transform
  input function is A[x][y][z], x in [0,nhx), y in [0,nhy), z in [0,nhz)
*/
{
  int x,y,z,nxy;
  fftcomplex *aux;

  nxy=max(nx,ny);
  allocarray(aux,nxy);

  loop (x,0,nx)
    loop (y,0,ny)
      Fourier(A[x][y],nz);

  loop (x,0,nx)
    loop (z,0,nz) {
      loop (y,0,ny) aux[y]=A[x][y][z];
      Fourier(aux,ny);
      loop (y,0,ny) A[x][y][z]=aux[y]; }

  loop (y,0,ny)
    loop (z,0,nz) {
      loop (x,0,nx) aux[x]=A[x][y][z];
      Fourier(aux,nx);
      loop (x,0,nx) A[x][y][z]=aux[x]; }

  free(aux);
}
