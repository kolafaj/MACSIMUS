#define StaNameLen 24 /* should be multiple of 8 */

void StaSet(double dt,int lag,int blklag,int noblk);
/*
  Set parameters for recording variables, 
  affects all subsequent first calls to StaAdd, cf. StaLag().
    dt = time between 2 pieces of calls to StaAdd
    lag = history length (in dt)
    blklag = history lengths for blocked data
    noblk = number of blockings, each by the factor of 2
*/

void StaAdd(const char *name, double x);
/*
  Record datum x under name name (max StaNameLen-1 chars)
*/

void StaFree(void);
/*
  Destroy all measurements and free the memory
*/

double StaPrint(const char *name,const char *how);
/*
  Print statistics of `name'. String `HOW' determines the mode:
    "/"    nothing printed
    NULL   only one line of mean+stderr is printed to out (statics=ONE)
    ""     compact table with c_t, tau, StDev printed to out (STA)
    "+"    as above with more decimal digits
    "-"    t:c_t table printed to stdout (TCF)
    "-@"   t:c_t table printed to file name.tcf (TCF)
    OTHER: t:c_t table printed to file HOW.name.tcf (TCF)
           (special characters in `name' are ignored - legacy mode)
    "-t"   in front of HOW = as above (TCF)
    "-T"   in front of HOW = as above + special chars kept except '/' -> '_'
    "-c"   in front of HOW = as -t except covariances are printed (COV)
    "-C"   in front of HOW = as above + special chars kept except '/' -> '_'

    REMOVED: "*" in front of HOW - is ERROR now!
*/

void StaPrintAll(const char *HOW);
/*
  Print statistics for all recorded variables, see above for `HOW`
*/

void StaSave(char *fn);
/*
  Save all measurements to file `fn'
*/

void StaKeySave(char *fn,int4 key);
/*
  As above and store one int4 value `key` at file end
  (to be honest, I do not know what is this good for)
*/

void StaLoad(char *fn);
/*
   Load all measurements saved by StaSave()
*/

extern int StaError;
/*
   Error indicator: set to 1 if `name' in the following functions is not found
*/

unsigned4 StaN(const char *name);
/*
  returns the number of measurements of `name'
*/

double StaMean(const char *name);
/*
  Returns the average of variable `name'
*/

double StaVar(const char *name);
/*
  Returns the variance of variable `name'
*/

double StaStdErr(const char *name);
/*
  Returns the estimate of the standard error of variable `name', 
  see StaErrMethod for the method
*/

double StaLag(const char *name,char key);
/*
  Returns selected propertes set by StaLag for `name' according to the key:
  'l': lag
  'k': blklag
  'n': noblk
  't': time = lag*DT
  'd': DT
  other: -9e9;
*/

void StaErrMethod(int m);
/*
  Set the method used to calculate the error estimate of the correlated
  time series, applies to subsequent call to StaStdErr, StaPrint, StaPrintAll
  0 = (default) The tables of stdev's contain sqrt(var*tau_n), where 
         tau_n=sum_i=1,n
      and c_i are time autocorrelation coefficients
      The final estimate is a heuristic average based on blocking and
      c_1,c_2 of the blocked data
  1 = The tables of stdev contain sqrt(var*tau), where tau for
      block-block uses the 'subaverage" formula for the 1st order 
      autoregressive process whereas higher correlations are extrapolated
      exponentially.  
      The final stdev uses the maximum available blocking + subaverage
      formula instead of 1+2c_1 (only a few blockings recommend)
  2 = As above.
      The final stdev uses the maximum available blocking + tau=1+2 c_1 +
      exponential extrapolation (1st order autoregressive process).
*/
