#define ARGLEN 102406

char *strend1(const char *str) /************************************ strend1 */
/* returns pointer to the last character (before 0) */
{
  char *c;

  for (c=(char*)str; *c; c++);

  return c-1;
}

char *dupstr(char *tok) /******************************************** dupstr */
{
  char *c=strdup(tok),*e;

  if (!c) Error("no heap");

  e=strend1(c);
  if (*e=='\n') *e=0;

  return c;
}

struct lett_s {
  struct lett_s *next;
  float x,y,r;
  char *s;
  char *cmd;
} *lett;

double number(char *c) /********************************************* number */
{
  char *end;
  double x=strtod(c,&end);

  if (end) {
    if (*end=='c') x*=72/2.54;
    else if (*end=='m') x*=72/25.4;
    else if (*end=='i') x*=72;
    else if (*end=='p') x*=72/72.27; }

  return x;
}

char *psdef="ps.def"; // or NAME.def
char *plotps="plot.ps"; // or NAME.ps
char *ploteps="plot.eps"; // or NAME.eps

void readpsdef(void) /******************************************** readpsdef */
/* caveat: memory leak */
{
  char line[ARGLEN];
  float rot=0;
  int istick=0;
  FILE *psdat=fopen(psdef,"rt");
  struct lett_s *l;
  static char *cmd;

  if (!psdat) return;

  /* remove old labeling */
  freelist(l,lett)

  while (fgets(line,ARGLEN,psdat)) {
    double x=number(line+2),y=0;
    char *c=strchr(line+3,' '),*cc=NULL,*parm2=NULL,*aa;
    int relx=line[2]=='+',rely=0;

    aa=strchr(line+1,' ');
    if (aa) aa++;      
    
    if (c) {
      c++;
      if (*c=='+') rely=1,c++;
      y=number(parm2=c);
      cc=strchr(c,' ');
      if (cc) cc++; }
      
    if (isdigit(line[0])) {
      /* redefine line */
      int col=atoi(line)%14-1;
      if (col<0) col=0;
      col=color[col];
      
      strcpy(strend1(line),",,,,,");
      
      cc=strchr(c=line+2,','); *cc++=0;
      my_style.line[col].rgb=dupstr(*c?c:"0 0 0");
      c=strchr(cc,','); *c++=0;
      my_style.line[col].dash=dupstr(cc);
      cc=strchr(c,','); *cc++=0;
      if (*c) my_style.line[col].thick=number(c);
      c=strchr(cc,','); *c++=0;
      if (*c) my_style.line[col].size=number(cc);
      cc=strchr(c,','); *cc++=0;
      my_style.line[col].symbol=dupstr(*c?c:"fullcircle"); }
    
    else switch (toupper(line[0])) {
        case 'B': {
          static struct my_style_line_s line[16] = {
            /*-*/    { "1 1 1", "", -1,"fullcircle",-1},
            /*14*/   { "0 0 0", "3 3 3 6", -1,"dotcircle",-1},
            /*12*/   { "0 0 0", "1 2 1 4 1 4 1 2", -1,"cross",-1},
            /*9*/    { "0 0 0", "1 2 1 2 1 4", -1,"fulldown",-1},
            /*13*/   { "0 0 0", "1 2 1 4 1 4", -1,"star",-1},
            /*10*/   { "0 0 0", "6 4 3 4", -1,"opendown",-1},
            /*8*/    { "0 0 0", "1 2 1 4", -1,"openup",-1},
            /*11*/   { "0 0 0", "1 4", -1,"plus",-1},
            /*-*/    { "0 0 0", "", -1,"fullcircle",-1},
            /*7*/    { "0 0 0", "6 3 1 3 1 3 1 3", -1,"fullup",-1},
            /*5*/    { "0 0 0", "3 4", -1,"fulldiamond",-1},
            /*3*/    { "0 0 0", "6 3 1 3", -1,"fullsquare",-1},
            /*6*/    { "0 0 0", "6 3 1 3 1 3", -1,"opendiamond",-1},
            /*4*/    { "0 0 0", "1 2", -1,"opensquare",-1},
            /*2*/    { "0 0 0", "6 4", -1,"opencircle",-1},
            /*1*/    { "0 0 0", "", -1,"fullcircle",-1} };
          copy(my_style.line,line,sizeof(line)); }
          break;

        case 'R':
          rot=x; break;

        case 'S':
          my_style.fontsize=x;
          if (parm2 && *parm2) {
            if (cc) { cc[-1]=0; my_style.enc=atoi(cc); }
            my_style.font=strdup(parm2); }
          break;

        case 'C':
          if (!aa) {
            cmd=NULL;
            break; }

          c=strchr(aa,'\n');
          if (c) *c=0;
          
          for (c=aa; *c; c++)
            if (isalpha(*c)) {
              cmd=strdup(aa); /* leak! */
              goto Break; }

          for (c=aa; *c; c++)
            if (isdigit(*c)) {
              alloc(cmd,strlen(aa)+12); /* leak! */
              strcpy(cmd,aa);
              strcat(cmd," setrgbcolor");
              goto Break; }
          cmd=NULL;
          Break:
          break;
          
        case 'L': {
          struct lett_s *l;
          static double lastx,lasty;

          allocone(l);
          l->next=lett;
          lett=l;
          l->x=x;
          if (relx) l->x=lastx+=x;
          else lastx=l->x=x;
          if (rely) {
            if (!y) y=my_style.fontsize*1.5;
            l->y=lasty-=y; }
          else
            lasty=l->y=y;
          lett->cmd=cmd;
          l->s=dupstr(cc); /* leak! */
          l->r=rot; }
          break;

        case 'X': 
          my_style.bx=x;
          if (c) my_style.bX=y;
          if (cc) my_style.xlabel=dupstr(cc);
          break;

        case 'Y': 
          my_style.by=x;
          if (c) my_style.bY=y;
          my_style.yrot=rot;
          rot=0;
          if (cc) my_style.ylabel=dupstr(cc); 
          break;

        case 'W': 
          /* negative = cm (legacy) */
          if (x<0) x*=-72/2.54;
          my_style.xsize=x;
          if (y<0) y*=-72/2.54;
          my_style.ysize=y;
          break;
        
        case 'M': my_style.mode=toupper(line[2]);
          if (my_style.mode=='P') {
            x=my_style.xsize;
            my_style.xsize=my_style.ysize;
            my_style.ysize=x; }
          break;

        case 'T': 
          my_style.thick=x;
          if (c) my_style.framethick=y;
          break;

        case 'F': 
          if (!istick++) my_style.nticks=0;
          if (my_style.nticks<5) {
            my_style.xticks[my_style.nticks]=my_style.yticks[my_style.nticks]=x;
            if (c) {
              my_style.yticks[my_style.nticks]=y;
              if (cc) my_style.tick[my_style.nticks]=number(cc); } }
          my_style.nticks++;
          break; 

        case '-': 
          my_style.minus=x; break;

        case '.': return;
      } 
  } 
}
