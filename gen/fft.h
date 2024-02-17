/* header file for fft.c, see also fourier.c */

#ifdef FLOAT
typedef float fftreal;
#else
typedef double fftreal;
#endif

typedef struct { fftreal re,im; } fftcomplex;

typedef struct Qstruct {
  struct Qstruct *next;
  int NofQ;
  fftcomplex Q[1]; } Qtype,*QPtr;

extern QPtr QHead;

/* returns pointer to the table of Q[j]==exp(i*j/n), j==0..n-1,
   creates new table if (necessary */
QPtr InitFour(int n);

/*
complex array a[0..m-1], m==abs(n), is replaced by the Fourier transform:

           m-1
    a[k] = SUM (a[j]*exp[i*j*k/n]), k==0..m-1
           j=1

NOTES:
* n<0 denotes back transform.
* Efficiency decreases if  n  contains large prime numbers,
  recommended values are  n == 2**k1 * 3**k2, where k2 is not too big. */
void Fourier(fftcomplex *a, int n);

/*
Fourier transform of 1 Even Real function a^ [0..nh] - not optimized !
Period of the function is n==2*nh.
'a^' may be a non-contiguous vector, increment is 'by' in units of
fftreal size (by==1 corresponds to a usual contiguous vector). */
void FourierER(fftreal *a, int nh, int by);

/* Fourier transform of 4 Even Real functions a^,b^,c^,d^ [0..nh].
  Optimized code, see comments of FourierER. */
void FourierER4(fftreal *a, fftreal *b, fftreal *c, fftreal *d, int nh, int by);

/*** version of ER3D optimized for x<-->y<-->z symmetric functions ***/
void Ff4(fftreal **p, int nh, int span);

void Fflush(fftreal **p, int nh, int span, int i);


/*
3D fourier transform of Even Real 3D function a[0..nh,0..nh,0..nh].
Period is 2*nh in each coordinate.
For a[x,y,z]=a[y,x,z]=a[y,z,x]...
NOTE map function: a[x,y,z] = x*span^2 + y*span + z, where span=nh+1
*/
void FourierER3D(fftreal *a,int nh);

/*
  3D general (complex) Fourier transform
  input function is A[x][y][z], x in [0,nhx), y in [0,nhy), z in [0,nhz)
*/
void Fourier3D(fftcomplex ***A,int nx,int ny,int nz);
