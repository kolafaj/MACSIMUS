/* #included twice: CALC=0,1 */

  if (*op & 128) switch (*op & 127) {
    /* unary operators (functions) */
    case '+': break;
    case '-': NR0=-NR0; break;
    case '\\': OP(nr,'\\'); NR0=sqrt(NR0); break;
    case 'u': OP(nr,'u'); NR0=cbrt(NR0); break;
    case 's': OP(nr,'f'); NR0=sin(NR0*_Id.degrad); break;
    case 'c': OP(nr,'f'); NR0=cos(NR0*_Id.degrad); break;
    case 't': OP(nr,'f'); NR0=tan(NR0*_Id.degrad); break;
    case 'S': OP(nr,'f'); NR0=asin(NR0)/_Id.degrad; break;
    case 'C': OP(nr,'f'); NR0=acos(NR0)/_Id.degrad; break;
    case 'T': OP(nr,'f'); NR0=atan(NR0)/_Id.degrad; break;
    case 'n': OP(nr,'f'); NR0= ~((unsigned)NR0); break; /* logical not */
    case 'b': {
      /* bit count */
      unsigned i=0,u=(unsigned)NR0;
      OP(nr,'f');
      while (u) { i+=u&1; u>>=1; } NR0=i; }
      break;
    case 'Q': OP(nr,'f'); NR0=NR0-(long int)NR0; break;
    case 'x': OP(nr,'f'); NR0=sinh(NR0); break;
    case 'y': OP(nr,'f'); NR0=cosh(NR0); break;
    case 'z': OP(nr,'f'); NR0=tanh(NR0); break;
    case 'X': OP(nr,'f'); NR0=asinh(NR0); break;
    case 'Y': OP(nr,'f'); NR0=acosh(NR0); break;
    case 'Z': OP(nr,'f'); NR0=atanh(NR0); break;
    case 'j': OP(nr,'f'); NR0=erf(NR0); break;
    case 'J': OP(nr,'f'); NR0=erfc(NR0); break;

    /* problems with tgamma,signgam; this is probably fool-proof: */
    case 'g': {
      int sg=1;
      double x=lgamma(NR0);
      OP(nr,'f');
      if (NR0<0) sg=(int)(NR0)&1?1:-1;
      NR0=exp(x)*sg;
      break; }
    case 'G': OP(nr,'f'); NR0=lgamma(NR0); break;
    case 'l': OP(nr,'f'); NR0=log(NR0); break;
    case 'd': OP(nr,'f'); NR0=log10(NR0); break;
    case 'w': OP(nr,'f'); NR0=log(NR0)/log(2.); break; /* lb()=log_2 */
    case 'e': OP(nr,'f'); NR0=exp(NR0); break;
    case 'B': NR0=fabs(NR0); break; /* abs */
      //    case 'u': OP(nr,'0'); NR0=NR0==0; break; /* null, nul */
    case 'I': OP(nr,'0'); if (NR0>0) NR0=1; else if (NR0<0) NR0=-1; break; /* sgn (used to be sign) */
    case 'H': OP(nr,'0'); if (NR0>0) NR0=1; else if (NR0<0) NR0=0; else NR0=0.5; break; /* Heaviside H() */
    case 'i': OP(nr,'f'); NR0=(int)(NR0); break;
    case 'F': OP(nr,'f'); NR0=floor(NR0); break;
    case 'A': /* float conversion: must cheat the optimization */
              { float F[2]; int i=1;
                F[1]=NR0;
                if (cos(2.3)>1) i=0;
                NR0=F[i]; }
              break;
#if CALC&1
    case 'v': { int i; loop (i,0,CALC_NUNITS) nr[0].pow[i]=0; } break;
    case 'V': NR0=1; break;
#endif
    case 'r': {

      OP(nr,'f');

      /* initialize */
      if (NR0<-1) {
        rndinit(7,(int)(-NR0));
        break; }
      
      if (!rndinit(-1,0)) {
        char *s=getenv("RNDSEED");
        int i=0;
        
        if (s) i=atoi(s);
        rndinit(7,i); }

      if (NR0==-1)
        /* [-1,1) */
        NR0=rndcos();
      else if (NR0==0)
        /* [0,1) */
        NR0=rnd();
      else if (NR0==1)
        /* normal Gaussian distribution */
	NR0=rndgauss();
      else
        /* integer in [0,NR0-1] */
	NR0=irnd(NR0);
      } break;
    /* factorial */
    case 'f':
      OP(nr,'f');
      if (NR0<0) goto reterror;
      if (NR0>1000) NR0=1000;
      { double f=1; while (NR0>0) f*=NR0,NR0-=1e0; NR0=f; }
      break;
    default: goto reterror; }

  else { nr--; switch (*op) {
    /* binary operators */
    case '&': OP(nr+1,'+'); NR0= (unsigned)NR0 & (unsigned)NR1; break;
    case '|': OP(nr+1,'+'); NR0= (unsigned)NR0 | (unsigned)NR1; break;
    case '$': OP(nr+1,'+'); NR0= (unsigned)NR0 ^ (unsigned)NR1; break;
    case '<': OP(nr+1,'<'); NR0= NR0 < NR1; break;
    case '>': OP(nr+1,'<'); NR0= NR0 > NR1; break;
    case 037: OP(nr+1,'<'); NR0= NR0 == NR1; break;
    case 036: OP(nr+1,'<'); NR0= NR0 >= NR1; break;
    case 035: OP(nr+1,'<'); NR0= NR0 <= NR1; break;
    case 034: OP(nr+1,'<'); NR0= NR0 != NR1; break;
    case '+': OP(nr+1,'+'); NR0+=NR1; break;
    case '-': OP(nr+1,'-'); NR0-=NR1; break;
    case '*':
      OP(nr+1,'*');
      if (NR0==0 || NR1==0) NR0=0; /* so that NaN*0=0 (2010-02-10) */
      else NR0*=NR1;
      break;
    case '/':
      OP(nr+1,'/');
      NR0/=NR1;
      break;
    case '@': /* because of ev.c - % is in format */
      //    case '`':  before 12/2023
    case '%':
      OP(nr+1,'f');
      NR0=NR0-((int)(NR0/NR1))*NR1;
      break;
    case '^':
      OP(nr+1,'^');
      NR0=pow(NR0,NR1); break;
    default:
      goto reterror; } }
