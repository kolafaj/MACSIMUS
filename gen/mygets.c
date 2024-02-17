/*** Replacement of the old good simple deprected gets()
 * reads one line from stdin (if not EOF)
 * getsbufsize>0:
   - stores max |getsbufsize| characters (incl. trailing 0)
   - prints a warning on overflow (extra characters are lost)
 * getsbufsize<0:
   - buffer overflow is error (exit(1))
 * returns:
   - NULL on EOF
   - the (optionally truncated) string read
***/
int getsbufsize=-1024;

char *mygets(char *s) /********************************************** mygets */
{
  int i,ch,ag=abs(getsbufsize);

  for (i=0; i<ag; i++) {
    ch=getchar();
    if (ch==EOF) { s[i]=0; return i?s:NULL; }
    if (ch=='\n') { s[i]=0; return s; }
    s[i]=ch; }

  s[ag-1]=0;
  if (getsbufsize>0)
    fprintf(stderr,"mygets() WARNING: string truncated because longer than %d\n",ag);
  else {
    fprintf(stderr,"mygets() ERROR: buffer %d overflow\n",ag);
    exit(1); }

  for (;;) {
    ch=getchar();
    if (ch==EOF || ch=='\n') return s; }
}

#define gets mygets
