/* files of variable record lengths - C version of VARFILE.PAS
                         (c) J.Kolafa 1990

NOTE: quite obsolete piece of software kept here for compatibility reasons
04/2016: #ifdef sdsalloc used instead of #ifdef SDS (sdsalloc must be macro)
         Reverse removed (see old+misc/varfile.h)
02/2009: extern variables merged to one structure for better (though still 
         cumbersome) handling of several open files
11/2001: VarError() and GetSize() made static to improve style

1 record requires (2*sizeof(int4)+size) bytes 
max length of a record is given by int4 (now 4-byte long)
record format: <int4 size><size bytes of data><int4 checksum>
size==0 denotes EOF 
NOTE VarGet no longer supported.
Since VarFile.size is known before the variable is read, you can use
  alloc(XXX,VarFile.size); VarRead(XXX,VarFile.size);  
*/

#include "ground.h"
#include "sds.h"
#include "varfile.h"

#define VarSum0 0x4A4B;

struct varfile_s VarFile;

/* local: */
static char VarOpenW;
static char VarErrStr[128];
static char *VarErrEnd=VarErrStr;

static void VarError(const char *msg,char *file,int line) /******** VarError */
{
  strcat(VarErrStr+strlen(VarErrStr),msg);
  extError(VarErrStr,file,line);
}
#undef Error
#define Error(X) VarError(X,__FILE__,__LINE__)

static void GetSize(void) /***************************************** GetSize */
{
  if (fread(&VarFile.size,sizeof(VarFile.size),1,VarFile.file)!=1) Error("size read");
  VarFile.eof = VarFile.size==0;
}

void VarOpen(char *fn, char *mode) /******************************** VarOpen */
{
  char m[3];

  sprintf(VarErrStr,"%sing \"%s\":",
          *mode=='r' ? "read" : *mode=='w' ? "writ" : "?-",
          fn);
  VarErrEnd=VarErrStr+strlen(VarErrStr);
  m[0] = *mode; m[1] = 'b'; m[2] = 0;
  if (!(VarFile.file=fopen(fn,m))) Error("open");
  VarOpenW = m[0]=='w';
  if ( ! (VarOpenW = m[0]=='w') ) GetSize();
}

void VarPut(void *v, int4 size) /************************************ VarPut */
{
  int4 i,s;
  unsigned char a;

  if (!size) Error("0 record to write");
  fwrite(&size,sizeof(size),1,VarFile.file);
  s = VarSum0;
  loop (i,0,size) {
    a = *((unsigned char*)v+i); 
    s += a; 
    putc(a,VarFile.file); };
  if (fwrite(&s,sizeof(s),1,VarFile.file)!=1) Error("write");
}

void VarReadN(char *name,void *v, int4 size) /********************* VarReadN */
{
  unsigned char a;
  int4 i,s;
  
  sprintf(VarErrEnd,"\nVarRead(%s,%d)",name,size);
  if (VarFile.eof) {
    Error("unexpected EOF"); return; }
  if (size != VarFile.size) {
    sprintf(VarErrEnd+strlen(VarErrEnd),", record length=%d",VarFile.size);
    Error(""); }
  s = VarSum0;
  loop (i,0,VarFile.size) {
    a = getc(VarFile.file); 
    *((unsigned char*)v+i) = a; 
    s += a; };
  if (fread(&i,sizeof(i),1,VarFile.file)!=1) Error("read");
  if (i!=s) Error("checksum");
  GetSize();
  *VarErrEnd=0;
}

#if 0
void *VarGet(void)
{
  void *v;
  
  if (v = (void *)malloc(VarFile.size)) { VarRead(v,VarFile.size); return v; }
  else Error("allocate");
  
  return v;
}
#endif

void VarClose(void) /********************************************** VarClose */
{
  int4 z=0;
  
  if (VarOpenW) {
    if (fwrite(&z,sizeof(z),1,VarFile.file)!=1) Error("write EOF"); }
  if (fclose(VarFile.file)) Error("close");
}

int4 KeyClose(int4 key) /****************************************** KeyClose */
{
  int4 z=0;
  
  if (VarOpenW) {
    if (fwrite(&z,sizeof(z),1,VarFile.file)!=1) Error("write EOF");
    if (fwrite(&key,sizeof(key),1,VarFile.file)!=1) Error("write key"); }
  else {
    if (key) if (fread(&key,sizeof(key),1,VarFile.file)!=1) Error("read key"); }
  if (fclose(VarFile.file)) Error("close");
  
  return key;
}

void VarPutSds(void *s) /***************************************** VarPutSds */
{
  VarPut((char*)s+sizeof(int4),sdssize(s)-sizeof(int4));
}

void *VarGetSds(void) /******************************************* VarGetSds */
{
  char *VarSds;
  
  sdsalloc(VarSds,VarFile.size+sizeof(int4));
  VarRead(VarSds+sizeof(int4),VarFile.size);
  
  return VarSds;
}

void VarReadSds(void *s) /*************************************** VarReadSds */
{
  VarRead((char*)s+sizeof(int4),sdssize(s)-sizeof(int4));
}
