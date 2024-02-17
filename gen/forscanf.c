/*
  fortran-like formatted input

  can be safely directly #included or used with forscanf.h as a module

  synopsis:

    int forscanf(FILE *file,const char *format,...);
    int sforscanf(char *line,const char *format,...);

  format field    corresponding parameter

    %9i %9d       pointer to integer
    %9li %9ld     pointer to long integer
    %9f           pointer to float
    %9lf          pointer to double
    %9a           array of char of length [9] (no \0 appended)
    %9s           array of char of length [9+1] (\0 appended, blanks removed)
    <space>       ignore space (fortran 1X)
    \n            ignore all data to end of input line

  the number after % is the length of the field irrespective of any
  blank characters (thus, spaces are included into strings!)

  max field length is 80

  tabulators are not and will not be supported (= they are treated as
  ordinary characters)

  number of scanned items is returned

  if \n is encountered before formats have been scanned, spaces are
  appended

  format %9a reads 9 characters

  format %9s reads 9 characters, removes leading and trailing spaces
  (' ', not tabs etc.) and appends \0; spaces inside the string are left
  (e.g., " A B " -> "A B" (neded by \0) if read by "%5s"

*/

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MYGETC0 { \
  c=fgetc(file); \
  if (c<0) { i=-1; goto done; } }
#define MYGETC { \
  if (eol) c='\n'; else c=fgetc(file); \
  if (c<0) { i=-1; goto done; } \
  if (c=='\n') { eol=1; c=' '; } }
int forscanf(FILE *file,const char *format,...)
{
#include "forscani.c"
}

#undef MYGETC
#undef MYGETC0
#define MYGETC0 c=line[pos++];
#define MYGETC { \
    if (eol) c='\n'; else { c=line[pos++]; if (!c) c='\n'; } \
  if (c=='\n') { eol=1; c=' '; } }
int sforscanf(char *line,const char *format,...)
{
  int pos=0;
#include "forscani.c"
}
