int yes(char *s)
{
  char line[256];
  printf("%s (y/n) ",s);
  gets(line);
  return line[0]=='y';
}

int GetInteger(char *s,int i)
{
  char line[256];
  printf("%s = %d ==> ",s,i);
  gets(line);
  if (line[0]) sscanf(line,"%d",&i);
  return i;
}

double GetReal(char *s,double r)
{
  char line[256];
  printf("%s = %g ==> ",s,r);
  gets(line);
  if (line[0]) sscanf(line,"%lf",&r);
  return r;
}
