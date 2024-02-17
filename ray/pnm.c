/* pnm read utilities - output is for depth=255 */

#define DITH 5
#define MPL 8
float Dither[DITH][DITH];

int xsize,ysize,depth;
double D[3];
int Depth[3]; /* passed depth */
double Gamma=1;
static int isgamma;
int dither[3];

FILE *pnm;
int pnmtype;
char comment[1024];
typedef unsigned char rgb_t[3];
rgb_t *rgb;
int rowsize,type;
unsigned char *row,*gray;
static int dbl;

int openpnm(char *fn)
/* suggested model returned (1: B/W, 2: G, 3: color) */
{
 char line[128];
 int i,j,k;
 int model;

 if (fn) pnm=fopen(fn,"rb"); else pnm=stdin;
 if (!pnm) Error(fn);

 isgamma=fabs(Gamma)!=1;

 loop (i,0,DITH) loop (j,0,DITH)
   Dither[i][j]=((i*DITH+j)*MPL%(DITH*DITH))/(float)(DITH*DITH);

 fgets(line,128,pnm);
 if (line[0]!='P') Error("bad header");
 type=atoi(line+1);
 if (type<1 || type>6) Error("unsupported file type, must be one of:\n\
text: P1 (bitmap), P2 (graymap), P3 (color pixmap)\n\
 bin: P4 (bitmap), P5 (graymap), P6 (color pixmap)");
 
 for (;;) {
   if (!fgets(line,128,pnm)) Error("unexpected EOF");
   if (line[0]!='#') break;
   else
     if (strlen(comment)+strlen(line)<1024) strcat(comment,line); }
 sscanf(line,"%d%d",&xsize,&ysize);
 if (strchr("\2\3\5\6",type)) {
   fgets(line,128,pnm);
   sscanf(line,"%d",&depth); }
 else depth=1;
 if (depth>255) {
   if (depth>65535) Error("depth>65535");
   fprintf(stderr,"WARNING: depth=%d reduced to depth=255\n",depth);
   dbl=1; }
 else 
   dbl=0;

 D[0]=D[1]=D[2]=255.99999/depth;
 loop (k,0,3) if (Depth[k]) {
   D[k]=255.99999/Depth[k];
   if (dbl) D[k]/=255.99999; }

 alloc(gray,xsize);
 alloc(rgb,xsize*3);
 
 switch (type) {
   case 4: case 1:
     model=1;
     rowsize=(xsize+7)/8;
     break;
   case 5: case 2:
     model=2;
     rowsize=xsize*(dbl+1);
     break; 
   case 6: case 3:
     model=3;
     rowsize=xsize*3*(dbl+1); 
     break;
   default: Error("unsupported P type"); }
 alloc(row,rowsize); 

 return model;
}

void readpnm(int inv,int j)
{
 int i,ii,k;
 static double w[3]={.299,.587,.114};

 if (type==1) loop (i,0,rowsize) {
   int k;
   row[i]=0;
   loop (k,0,8) {
     if (i*8+k>=xsize) ii=0;
     else
       if (1!=fscanf(pnm,"%d",&ii)) Error("unexpected EOF");
     row[i]=(row[i]<<1)|ii; } }
 else if (type<4) loop (i,0,rowsize) {
   if (1!=fscanf(pnm,"%d",&ii)) Error("unexpected EOF");
   row[i]=ii; }
 else
   if (1!=fread(row,rowsize,1,pnm)) Error("unexpected EOF");

 switch (type) {
   case 4: case 1:
     loop (i,0,xsize) {
       ii=!!(row[i/8]&(128>>(i&7)));
       if (inv) switch (inv) {
         case 1: case 4: case -5: case -6: ii=!ii; break;
         case -1: inv=ii=!ii; }
       rgb[i][0]=rgb[i][1]=rgb[i][2]=gray[i]=ii*255; }
     break;
   case 5: case 2:
     loop (i,0,xsize) {
       double dii;

       if (inv) {
	 if (inv==-1) inv=((short unsigned*)row)[i]<248;
	 if (inv==4 || inv==6 || inv==-5) inv=0; }

       if (dbl) 
         dii=((short unsigned*)row)[i]*D[0]+(dither[0]*Dither[i%DITH][j%DITH]);
       else 
         dii=row[i]*D[0]+(dither[0]*Dither[i%DITH][j%DITH]);
       Min(dii,255.)
       if (inv) dii=255.-dii;
       if (isgamma)
         if (Gamma>0) dii=255*pow(dii/255.,Gamma)+0.499;
         else dii=255*(1-pow(1-dii/255.,-Gamma))+0.499;
       rgb[i][0]=rgb[i][1]=rgb[i][2]=gray[i]=(int)(dii+0.5); }
     break;
   case 6: case 3:
     loop (i,0,xsize) {
       double dd[3],q=0;

       if (inv) {
	 if (inv==-1) inv=((short unsigned*)row)[3*i]*(int)30+((short unsigned*)row)[3*i+1]*59+((short unsigned*)row)[3*i+2]*11<24800;
	 if (inv==4 || inv==5 || inv==-6) inv=0; }

       loop (k,0,3) {
         if (dbl) 
           dd[k]=((short unsigned*)row)[3*i+k]*D[k]+(int)(dither[k]*Dither[i%DITH][j%DITH]);
         else
           dd[k]=row[3*i+k]*D[k]+(int)(dither[k]*Dither[i%DITH][j%DITH]);
         Min(dd[k],255.) 
         if (inv) dd[k]=255.-dd[k];
         q+=w[k]*dd[k]; }

       if (isgamma)
         if (Gamma>0) {
           q=pow(q/255.,Gamma-1);
           loop (k,0,3) dd[k]*=q; }
         else {
           q=pow(1-q/255.,Gamma)/q;
           loop (k,0,3) dd[k]*=q; }

       q=0;
       loop (k,0,3) {
         q+=dd[k]*w[k];
         ii=(int)(dd[k]+0.5);
         Min(ii,255)
         rgb[i][k]=ii; }
       ii=(int)(q+0.5);
       Min(ii,255)
       gray[i]=ii; }
   }
}

