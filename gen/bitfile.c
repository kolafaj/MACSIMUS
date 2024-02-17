#include "ground.h"
#include "bitfile.h"

struct bitfile_s bitfile={1};

int ceillog2(unsigned4 n) /*************************************** ceillog2 */
/*
  returns # of bits necessary to hold number n
*/
{
int i=1;
while (n>>=1) i++;
return i;
}

static byte getbyte(void) /**************************************** getbyte */
/*
  gets one byte from file bitfile.f - to be eventually buffered
*/
{
byte B;
bitfile.err=fread(&B,1,1,bitfile.f)!=1;
bitfile.bytecount++;
return B;
}

static void putbyte(byte B) /************************************** putbyte */
/*
  puts one byte to file bitfile.f - to be eventually buffered
*/
{
bitfile.err=fwrite(&B,1,1,bitfile.f)!=1;
bitfile.bytecount++;
}

void assignbit(FILE *f,char mode) /****************************** assignbit */
{
bitfile.f=f;
if (!strchr("rw",mode)) Error("bitfile open mode");
bitfile.err=0;
bitfile.mode=mode;
bitfile.bytecount=0;
bitfile.checksum=0;
if (bitfile.pos) Error("bitfile re-assigned and not flushed");
}

unsigned4 getbit(int b) /******************************************* getbit */
/*
  gets b bits from file bitfile.f
*/
{
unsigned4 val;
while (b>bitfile.pos) {
  bitfile.buf = (bitfile.buf<<8) | getbyte();
  bitfile.pos += 8; }
bitfile.pos -= b;
val=(bitfile.buf>>bitfile.pos) & ((1L<<b)-1);
bitfile.checksum+=val;
return val;
}


void putbit(unsigned4 val,int b) /********************************** putbit */
/*
  puts b bits of val to file bitfile.f, if check then checks whether
  val fits into b bits
  max b is 24 bits
*/
{
if (bitfile.check) {
  if (ceillog2(val)>b)
    ERROR(("fix overflow: %lu has more than %d bits",(long unsigned)val,b))
  if (b>24) ERROR(("fix overflow (%d is more than 24 bits)",b)) }
bitfile.checksum+=val;
while (bitfile.pos>=8) {
  bitfile.pos -= 8;
  putbyte((bitfile.buf>>bitfile.pos)&0xff); }
bitfile.buf = (bitfile.buf<<b) | val;
bitfile.pos += b;
}

void flushbit(void) /********************************************** flushbit */
/*
  mode='w': flushes buf into bitfile.f: next write will start from new byte
  mode='r': moves pointer to byte boundary (to be used on places 
            file was written with flushbit('w') )
*/
{
if (bitfile.mode=='w') while (bitfile.pos>0) {
  bitfile.pos -= 8;
  putbyte(
    (bitfile.pos>=0 ? bitfile.buf>>bitfile.pos : bitfile.buf<<-bitfile.pos)
    &0xff ); }
bitfile.pos=0;
}

