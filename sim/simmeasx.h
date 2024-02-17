
// void diffSF(double dt,double dtcp,int no,int n0,int nspec);

void initSF(void);
void calculateSF(void);
void printfSF(void);

void initdiff(double dtplb);
void calculatediff(int n,int no);
void printdiff(void);

#ifdef XSECTION
double mXsection(ToIntPtr A);
double cXsection(ToIntPtr A);
#  ifdef CLUSTERS
double clXsection(ToIntPtr A);
#  endif /*# CLUSTERS */
#endif /*# XSECTION */

#ifdef CLUSTERS
/* user interface, see the manual for details */
typedef struct cl_s {
  int mode; /* bits: 1=on, 2=clusters, 4=configurations, 8=bonddynamics */
  int format; /* bits: 1=full list, 2=alt list, 
                       4=cluster sizes, 8=alt sizes, 
                       16=prepend time, 32=dump clusters >= cl.mincluster */
  int maxn; /* If given, cluster counts up to size maxn are printed to the
               prt-file for each frame analyzed in the format defined by
               variable cl.format.
               Cluster topology is not distinguished here, only sizes matter. */
  int maxcluster; /* larger clusters are not topologically distinguished */
  int mincluster; /* ignore all smaller in detailed output (format&32) */
} cl_t;
extern cl_t cl;

void readclusterdef(char *fn);
void analyzeclusters(int frame);
void printclusters(void);
#endif /*# CLUSTERS */

#ifdef BJERRUM
void bjerrum_read_oo(void);
void bjerrum_topology(ToIntPtr A);
void bjerrum_gauss(ToIntPtr A,double from,double to,double q,double eps);
#endif

void measuredrifts(void);
#ifdef SLAB
void measure1drift(int n);
#endif
