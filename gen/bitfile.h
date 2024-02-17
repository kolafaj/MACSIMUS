typedef unsigned char byte;

extern struct bitfile_s {
  int check;           /* whether putbit should check fixed overflow */
  int pos;             /* bit position in buffer buf */
  unsigned4 buf;       /* bit i/o buffer */
  FILE *f;             /* file of bits */
  unsigned4 bytecount; /* # of bytes transfered */
  unsigned4 checksum;  /* checksum */
  int err;             /* set if r/w error on f */
  char mode;           /* 'r'/'w' */
  } bitfile;

void assignbit(FILE *f,char mode); /* assign file f for bit input (mode='r')
				      or output ('r') */
unsigned4 getbit(int b);           /* read b bits */
void putbit(unsigned4 val,int b);  /* write b bits (max 24) */
void flushbit(void);               /* flush bit buffer to file f 
                                      (does NOT flush file to disk!) */

int ceillog2(unsigned4 n);         /* # of bits necessary to hold number n */
