/*
04/2016: #includable module sds.h instead of #define SDS for alloc.h

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Self-defined structures (SDS) contain their size as int at the beginning
% (the size of this int is included into the count)
%
% If A and B are pointers to SDS, then
%   sdsalloc(A,<size in bytes>)  allocates SDS A.  The size is also written to
%                                the first int
%   sdsralloc(A,<size in bytes>) as above with ralloc
%   sdsalloczero,sdsralloczero : as above plus fills the memory by zeroes
%   sdscopy(A,B)  copies *B into *A.  Length (but not overlap...) are checked
%   sdszero(A) fills the contents of *A by zeroes (excl. the size)
%
% example:
%   typedef struct { int size; double a[10]; } xtype;
%   xtype X={sizeof(X)}, *A;
%   sdszero(&X);
%   sdsalloc(A,sizeof(xtype));
%   sdscopy(A,&X);
%
% to be used as include file
% required: 
%   macro Error(char *msg) and "ground.c" (originally "alloc.c")
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#  define sdscopy(X,Y) { \
  if (*(int*)(X) < *(int*)(Y)) ERROR(("SDS:lengths")) \
  else copy((int*)(X)+1,(int*)(Y)+1,*(int*)(Y)-sizeof(int)); }
/*.....  else if (abs((char*)(X)-(char*)(Y)) < *(int*)(X)) Error(sdserrovl); \*/

#  define sdszero(X)  { \
  if (*(int*)(X) < sizeof(int)) ERROR(("SDS:bad length")) \
  else memset((int*)(X)+1,0,*(int*)(X)-sizeof(int)); }

#  define sdssize(X) (*(int*)(X))

#  define sdsalloc(X,Y) { alloc((X),Y); *(int*)(X)=Y; }
#  define sdsalloczero(X,Y) { alloczero((X),Y); *(int*)(X)=Y; }

#  define sdsralloc(X,Y) { ralloc((X),Y); *(int*)(X)=Y; }
#  define sdsralloczero(X,Y) { ralloczero((X),Y); *(int*)(X)=Y; }
