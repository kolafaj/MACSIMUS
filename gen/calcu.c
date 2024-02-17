/* unit order: m,kg,s,K,mol,A */
struct unit_s Nr[NRSTACKLEN],*nr,retzero;

/* MACSIMUS program units (p.u.) */
static double PU[CALC_NUNITS] = {
/*  m         kg          s    K       mol                A */
  1e-10, 1.380649e-27,  1e-12, 1, 1/6.02214076e23, 3.919412183482151e-10 };

static struct unit_s unitlist[] = {
  {"m",1,{1}},
  {"kg",1,{0,1}}, /* bug: kkg allowed */
  {"s",1,{0,0,1}},
  {"K",1,{0,0,0,1}},
  {"mol",1,{0,0,0,0,1}},
  {"A",1,{0,0,0,0,0,1}},
  {"g",1e-3,{0,1}}, /* to parse mg etc. */
  {"ct",0.2e-3,{0,1}}, /* carat */
  {"mcg",1e-6,{0,1}}, /* ug */
  //  {"cd",1,{0,0,0,0,0,0,1}},
  //  {"sr",1,{0,0,0,0,0,0,0,1}},

  //  {"%",1e-2,{0}}, // too difficult to separate %
  {"ppm",1e-6,{0}},
  {"ppb",1e-9,{0}},

  {"B",8,{0}},

  {"au",149597870700,{1}}, /* astronomical unit, before u */
  {"pc", 30.856775814913673e15, {1}}, /* parsec */

  {"J",1,{2,1,-2}},
  {"W",1,{2,1,-3}},
  {"Wh",3600,{2,1,-2}},
  //  {"Ry",2.179872324923e-18,{2,1,-2}}, /* Rydberg (energy) */
  //  {"Rinf",10973731.568508,{-1}}, /* Rydberg (inverse length) */
  {"N",1,{1,1,-2}},
  {"Pa",1,{-1,1,-2}},
  {"P",0.1,{-1,1,-1}}, // Poise = g.cm-1.s-1
  {"C",1,{0,0,1,0,0,1}},
  {"V",1,{2,1,-3,0,0,-1}}, // volt = m2·kg·s−3·A−1
  {"S",1,{-2,-1,3,0,0,2}}, // siemens = kg−1·m−2·s3·A2
  {"Mho",1,{-2,-1,3,0,0,2}}, // siemens = kg−1·m−2·s3·A2
  {"Ohm",1,{2,1,-3,0,0,-2}}, // Ohm = kg·m2·s-3·A-2
  //  {"Ω",1,{2,1,-3,0,0,-2}}, // Ω = kg−1·m−2·s3·A2
  {"F",1,{-2,-1,4,0,0,2}}, // Farad = 1 s4·A2·m−2·kg−1
  {"H",1,{2,1,-2,0,0,-2}}, // 1 Henry = kg⋅m2⋅s-2⋅A-2
  {"T",1,{0,1,-2,0,0,-1}}, // 1 Tesla = kg·s−2·A−1
  {"G",1e-4,{0,1,-2,0,0,-1}}, // 1 Gauss = 1e-4 Tesla
  {"Wb",1,{2,1,-2,0,0,-1}}, // 1 Weber = 1 kg·m2·s-2·A-1
  //  {"D",1e-21/299792458,{1,0,1,0,0,1}}, Debye old SI
  // sqrt(2e-49[N m4]*pi*e^2/(alpha*h*c)) new SI
  {"D",3.3356409510722084e-30,{1,0,1,0,0,1}},
  {"ar",1e2,{2}}, /* ar - not a which collides with a=annus */
  {"ha",1e4,{2}}, /* hectare (accepted as non-SI unit by CIPM and EEC), before a */
  {"min",60,{0,0,1}},
  {"h",3600,{0,0,1}},
  {"d",3600*24,{0,0,1}},
  {"Ci",3.7e10,{0,0,-1}},
  {"Hz",1,{0,0,-1}},
  {"RPM",1./60,{0,0,-1}},
  {"a",  31556925.445,{0,0,1}}, /* annus = tropical year 2000 */
  {"a_T",31556925.216,{0,0,1}}, /* annus = mean tropical year */
  {"a_J",31557600,{0,0,1}}, /* annus = Julian year */
  {"a_G",31556952,{0,0,1}}, /* annus = Gregorian year */
  {"Bq",1,{0,0,-1}},
  {"L",1e-3,{3}},
  {"gal",3.785411784e-3,{3}}, /* 231 in^3 */
  {"USgal",3.785411784e-3,{3}},
  {"UKgal",4.54609e-3,{3}},
  {"M",1e3,{-3,0,0,0,1}}, /* mol/L */
  {"Da",1.6605390671738466e-27,{0,1,0,0,0}}, /* "mass-based-interpretation" */
  {"u",1.6605390671738466e-27,{0,1,0,0,0}}, /* "mass-based-interpretation" */
  {"kat",1,{0,0,-1,0,1}}, /* katal */
  {"AA",1e-10,{1,0,0,0,0,0}},
  {"eV",1.602176634e-19,{2,1,-2}},
  {"cal",4.184,{2,1,-2}},
  {"calIT",4.1868,{2,1,-2}},
  {"calIUNS",4.182,{2,1,-2}},
  {"in",0.0254,{1}},
  {"ft",12*0.0254,{1}},
  {"yd",36*0.0254,{1}},
  {"mi",1760*36*0.0254,{1}},
  {"nmi",1852,{1}},
  {"BTU",1055.056,{2,1,-2}},
  {"HP",745.69987158227,{2,1,-3}}, /* mechanical, imperial */
  //  {"HPe",746,{2,1,-3}}, /* electric */
  {"PS",735.49875,{2,1,-3}}, /* metric */
  {"lb",0.45359237,{0,1}},
  {"PSI",6894.757293168,{-1,1,-2}},
  {"bar",1e5,{-1,1,-2}},
  {"atm",101325.,{-1,1,-2}},
  {"at",98066.5,{-1,1,-2}},
  {"torr",101325./760,{-1,1,-2}},
  {"foe",1e44,{2,1,-2}}, /* 10 to fifty-one erg */

  {"t",1e3,{0,1}}, /* metric ton - must be after "ft", "ct" */
  {"erg",1e-7,{2,1,-2}},
  {"dyn",1e-5,{1,1,-2}},
  {"b",1e-28,{2}}, /* barn */
  //  {"a",1e2,{2}}, /* NOT accepted by SI, collision with a=annus */

  {NULL,0,{0}}};

static struct pref_s {
  char *name;
  double prefix;
} prefixes[] = {
  {"",1},
  {"q",1e-30}, // quecto
  {"r",1e-27}, // ronto
  {"y",1e-24}, // yocto
  {"z",1e-21}, // zepto
  {"a",1e-18}, // atto
  {"f",1e-15}, // femto
  {"p",1e-12}, // pico
  {"n",1e-9},  // nano
  {"u",1e-6},  // micro
  //  {"μ",1e-6},
  {"m",1e-3},  // milli
  {"c",1e-2},  // centi
  {"d",1e-1},  // deci
  {"da",10},   // deca
  {"h",100},   // hecto
  {"k",1e3},   // kilo
  {"ki",1024.},
  {"M",1e6},   // mega
  {"Mi",1024.*1024.},
  {"G",1e9},   // giga
  {"Gi",1024.*1024.*1024.},
  {"T",1e12},  // tera
  {"Ti",1024.*1024.*1024.*1024.},
  {"P",1e15},  // peta
  {"E",1e18},  // exa
  {"Z",1e21},  // zetta
  {"Y",1e24},  // yotta
  {"R",1e27},  // ronna
  {"Q",1e30},  // quetta
  {NULL,0}};

void towards(char *u);
char *to_remove; /* postponed removal of an unknown preferred unit */

static
void parseunit(struct unit_s *nr,char *b,char **end,int df) /***** parseunit */
/*
   parses "unit]", return value in nr->pow[]
   returns pos after ] in *end
   df=0 : does not store the unit into lastpow
   df=1 : stores the unit into lastpow (to be used for [])
   examples:
     [ms]=1e-3 s
     [m s]=m*s
     [m1/2]=[m0.5]=m^0.5
     [kg/m s]=[kg/m/s]=kg/(m*s)
     [us]=1e-6 s
     [m.s] = [m s]
(    [m.s]=m^0*s (because "." expands to "0.0", which is zero) - FIXED )
*/
{
  int i,j,k,sg=1;
  char *a,*e;
  double prefix=1;
  int pu;
  static double lastpow[CALC_NUNITS];
  static double lastprefix=1;

  /* MACSIMUS program unit */
  if (*b=='-') pu = -1;
  else if (*b=='+') pu = 1;
  else pu = 0;
  b += abs(pu);

  /* [] = previous unit */
  if (*b==']' || *b==0) {
    *end=b+(*b==']');
    loop (k,0,CALC_NUNITS) nr->pow[k]=lastpow[k];
    nr->val*=lastprefix;
    return; }

  loop (k,0,CALC_NUNITS) nr->pow[k]=0;

  /* trailing ] may be missing */
  e=strchr(b,']');
  if (!e) e=b+strlen(b)-1;

  while (*b && *b!=']') {
    double ex=1; /* exponent of unit */

    while (*b && *b<=' ') b++;
    /* test divide by unit (not part of fractional exponent)
       [m/s/kg]=[m s-1 kg-1] */
    if (*b=='/') sg=-1,b++;
    /* this version would interpret [m/s/kg]=[m s-1 kg]:
       if (*b=='/') sg=-sg,b++; */
    while (*b && *b<=' ') b++;

    if (!isalpha(*b)) {
      fprintf(stderr,"calcu: syntax [%s\n",b);
      return; }

    for (a=b; isalpha(*a) || *a=='_'; a++);
    /* unit (w. prefix) = b to a-1 incl. */

    /* optional exponent expected at a */
    ex=strtod(a,&e);
    if (e==a)
      /* no exponent -> 1 */
      ex=1;
    else {
      /* exponent found */
      while (*e && *e<=' ') e++;
      if (*e=='/') {
        /* fractional exponent ? */
        char *ee;
        double den;

        den=strtod(e+1,&ee);
        if (e+1!=ee) {
          /* yes, fractional exponent */
          ex/=den;
          e=ee; } } }

    for (i=0; unitlist[i].name; i++) {
      char *unitbeg=a-strlen(unitlist[i].name); /* if found */

      if (!memcmp(unitbeg,unitlist[i].name,strlen(unitlist[i].name))) {

        j=0; /* no prefix -> 1 */
        if (unitbeg==b) goto doprefix;
        if (unitbeg>b)
          for (j=1; prefixes[j].name; j++)
            if (!memcmp(b,prefixes[j].name,unitbeg-b)) goto doprefix;
        /* prefix not found, trying name (e.g., atm: is not at+m) */
        continue;
      doprefix:
        loop (k,0,CALC_NUNITS) nr->pow[k]+=unitlist[i].pow[k]*ex*sg;
        if (pu>=0) prefix*=pow(unitlist[i].val*prefixes[j].prefix,ex*sg);
        goto found; } }

    fprintf(stderr,"calcu: unknown unit [%s (ignored)\n",b);

    to_remove=malloc(strlen(b)+4);
    sprintf(to_remove," ~%s",b);
    
    /* not found - ignore for now */
  found:;
    /* NB: syntax with a dot [J.s] solved here: */
    while (*e && (*e<=' ' || *e=='.')) e++;
    b=e; }

  if (!*e) *end=e; /* missing ] */
  else *end = e+1;

  if (pu<0)
    /* MACSIMUS units */
    loop (k,0,CALC_NUNITS) prefix*=pow(PU[k],nr->pow[k]);
  else if (pu>0)
    loop (k,0,CALC_NUNITS) {
      prefix/=pow(PU[k],nr->pow[k]);
      nr->pow[k]=0; }

  nr->val*=prefix;

  /* store last unit */
  if (df) {
    loop (k,0,CALC_NUNITS) lastpow[k]=nr->pow[k];
    lastprefix=prefix; }

  return;
}

char *unitserr;

void opunits(struct unit_s *a,struct unit_s *b,char assignopb) /**** opunits */
/*
   combines units of two operands, checks if permitted
   a=1st operand and result, b=2nd operand
   E.g.: [m]+[kg] is invalid, [m]*[kg]=[m kg], [m]=[m] is [1],...
*/
{
  int k;
  unitserr=NULL;

  loop (k,0,CALC_NUNITS) switch (assignopb) {
    case 'p':
      /* parm[]: check */
      if (a->pow[k]!=b->pow[k]) unitserr="units of parameters do not match\n";
      break;
    case '+':
    case '-':
      /* sum, etc.: check if units match */
      if (a->pow[k]!=b->pow[k]) unitserr="cannot add,subtract different units\n";
      break;
    case '<':
    case '>':
      /* relation operators: check + set dimensionless */
      if (a->pow[k]!=b->pow[k]) unitserr="cannot compare different units\n";
      a->pow[k]=0;
      break;
    case 'f':
      /* check no unit */
      if (a->pow[k]!=0) unitserr="function argument is not dimensionless\n";
      break;
    case 'i':
      /* check no unit */
      if (a->pow[k]!=0) unitserr="BUG: iter,repeat,contfrac require dimensionless args\n";
      break;
    case '*': a->pow[k]+=b->pow[k]; break;
    case '/': a->pow[k]-=b->pow[k]; break;
    case '\\': a->pow[k]/=2; break; // square root (sqrt)
    case 'u': a->pow[k]/=3; break; // cube root (cbrt)
    case '^':
      /* power */
      if (b->pow[k]!=0) unitserr="exponent in power is not dimensionless\n";
      a->pow[k]*=b->val;
      break;
    case '0':
      /* set zero */
      a->pow[k]=0;
      break;
    default: fprintf(stderr,"internal: bad opunits key\n");
  }

  if (unitserr) fputs(unitserr,stderr);
}

/* user-friendly units - if found, this unit will be used */
struct user_s {
  struct user_s *next;
  char *unit;
} *userhead;

void removeunit(struct user_s *u) /****************************** removeunit */
{
  struct user_s *uu;

  if (!u) return;

  if (u==userhead) {
    userhead=userhead->next;
    free(u);
    return; }
  else
    looplist (uu,userhead) if (uu->next==u) {
      uu->next=u->next;
      free(u);
      return; }
}

void towards(char *u) /********************************************* towards */
/* command to, u points to the argument
  register the preferred form of unit: to J/s, to [J/s]
  remove it: to ~J/s, to ~[J/s]
  remove all: to ~
 */
{
  struct user_s *user;
  int rm=0;
  char *ket;
  //  fprintf(stderr,"/// %s///\n",u);
  while (*u==' ' || *u=='\t') u++;

  if (*u==0) return;

  if (*u=='~') rm=1,u++;
  while (*u==' ' || *u=='\t') u++;

  if (*u==0) { /* remove all */
    struct user_s *uu,*next;

    for (uu=userhead; uu; uu=next) {
      next=uu->next;
      free(uu); }
    userhead=NULL;

    return; }

  if (*u=='[') u++;

  ket=strchr(u,']');
  if (ket) *ket=0;

  if (rm) {
    looplist (user,userhead) if (!strcmp(user->unit,u)) {
        removeunit(user);
        return; }
    /* no error message if not found */
    return; }

  user=malloc(sizeof(*user));
  if (!user) { fprintf(stderr,"no heap\n"); exit(1); }
  user->next=userhead;
  userhead=user;
  user->unit=strdup(u);
}

double unitfactor=1;

char *unit(struct unit_s *res) /*************************************** unit */
/* return unit as a string */
{
  static char prtunit[128];
  char *p=prtunit+1;
  int k;
  int nz=0;
  struct unit_s nr;
  struct user_s *u;

  looplist (u,userhead) {
    char *end;

    //    fprintf(stderr,"u->unit=%p %s\n",u,u->unit);
    prtunit[0]='['; strcpy(p,u->unit); strcat(p,"]");
    nr.val=1;
    parseunit(&nr,p,&end,0);
    //    fprintf(stderr,"parseunit %s)\n",p);
    loop (k,0,CALC_NUNITS) if (res->pow[k]!=nr.pow[k]) goto notfound;
    unitfactor=nr.val;
    return prtunit;
  notfound: ; }

  unitfactor=1;
  *p=0;
  loop (k,0,CALC_NUNITS) if (res->pow[k]) {
    nz++;
    if (res->pow[k]==1)
      p+=sprintf(p,"%s ",unitlist[k].name);
    else
      p+=sprintf(p,"%s%g ",unitlist[k].name,res->pow[k]); }
  if (nz) {
    prtunit[0]='[';
    p[-1]=']'; }
  else
    prtunit[0]=0;

  return prtunit;
}

double unitf(struct unit_s *res) /************************************ unitf */
{
  unit(res);
  return unitfactor;
}
