char *posA,*negA;

char *optA(char *c) /************************************************* optA */
{
c+=2;
while (strchr("-0123456789",*c)) c++;
return c;
}

int groupcharge(int c) /**************************************** groupcharge */
/* in % */
{
double gc=site[c].charge;
int i;

loop (i,0,nc) {
  int b0=bond[i][0],b1=bond[i][1];
  if (b0!=b1) if (c==b0 || c==b1) gc+=site[c==b0?b1:b0].charge; }

return (int)(gc*100.+9000.5)-9000;
}

