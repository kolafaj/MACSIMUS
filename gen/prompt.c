/* 
   line editing with prompt, directly #include, link with -lncurses
   used in: c/cur.c c/ev.C
   possible setup (see c/ev.C):
     initscr(); cbreak();
     noecho(); nl();
     scrollok(stdscr,TRUE);
     intrflush(stdscr, FALSE);
   stop by:
     endwin()
   uses:
     getch() = read 1 char
*/

#include <curses.h>

#define PROMPT_N 1024 // line length
#define PROMPT_H 256 // history length

struct {
  char prompt[8]; 
  int n; /* line length */
  int pos; /* cursor position */
  int point; /* point from which the line is shown (if scrolled) */
  int back;
  char hist[PROMPT_H][PROMPT_N];
} _pr = { "", 79, 0, 0, 1 };

#ifdef DEBUG
FILE *debug;
#endif

void prompt_write(char *fn) /********************************** prompt_write */
{
  FILE *f=fopen(fn,"wt");
  int i;

  for (i=PROMPT_H-1; i>=0; i--) if (strlen(_pr.hist[i])>2)
    if (memcmp(_pr.hist[i]+2,"write ",6)) if (strcmp(_pr.hist[i]+2,"write"))
      fprintf(f,"%s\n",_pr.hist[i]+2);
  fclose(f);
}

void prompt_showline(char *line) /************************** prompt_showline */
{
  int pos0=_pr.pos,i,l=strlen(line);

  if (_pr.pos<=_pr.point) _pr.point=_pr.pos-_pr.n/6;
  if (_pr.pos>=_pr.point+_pr.n-1) _pr.point=_pr.pos-_pr.n*5/6;
  if (_pr.point<0) _pr.point=0;
  addch(13);
  for (i=_pr.point; i<_pr.point+_pr.n; i++)
    addch(i==_pr.point&&_pr.point ? '~' :
          i==_pr.point+_pr.n-1&&i<l-1 ? '~' :
          i<l ? line[i] : ' ');

  if (_pr.back) {
    _pr.pos=_pr.point+_pr.n;
    while (_pr.pos>pos0) addch(8),_pr.pos--;
    refresh(); }
}
    
int prompt_reshow(char *line) /******************************* prompt_reshow */
{
  if (_pr.pos<_pr.point || _pr.pos>=_pr.point+_pr.n) {
    prompt_showline(line);
    return 1; }
  return 0;
}

char *prompt_getline(char *line) /*************************** prompt_getline */
{
  int l=strlen(_pr.prompt),h=-1;
  char *ch,c;
  static int insmode=1,lastc;

  memset(line,0,PROMPT_N);
  strcpy(line,_pr.prompt);
  _pr.pos=strlen(line);
  prompt_showline(line);

  do {
    lastc=c;
    c=getch();

#ifdef DEBUG
    fprintf(debug,"%d\n",c); fflush(debug);
#endif

    if (c>=' ' && c<=126) { 
      if (insmode && _pr.pos<strlen(line)) {
        char a=line[_pr.pos],b;
        
        for (ch=line+_pr.pos; *ch; ch++) {
          b=ch[1]; ch[1]=a; a=b; } 
        prompt_showline(line); }
      addch(line[_pr.pos++]=c); 
      prompt_reshow(line); }

    if (c==27)
      switch (getch()) {
        /* esc-O = xterm-style, esc-[ = gnome-term style */
        case 'O':
        case '[': switch (getch()) {
          case '2': 
            /* Insert */
            if (getch()=='~') insmode=!insmode;
            break;
          case '3':
            /* Delete */
            if (getch()=='~') {
              for (ch=line+_pr.pos; *ch; ch++) *ch=ch[1];
              *ch=ch[l];
              prompt_showline(line); }
            break;
          case 'H': /* HOME */
            _pr.point=0, _pr.pos=l;
            prompt_showline(line);
            break;
          case 'F':  /* END */
            _pr.pos=strlen(line);
            prompt_showline(line);
            break;
          case 'A': /* up */
            if (h<PROMPT_H-1 && _pr.hist[h+1][0]!='?') h++; 
            goto B;
          case 'B': /* down */
            if (h) h--;
           B:
            if (h>=0 && h<PROMPT_H) {
	      memset(line,0,PROMPT_N);
	      strcpy(line,_pr.hist[h]); }
            _pr.pos=strlen(line); _pr.point=0;
            prompt_showline(line);
            break;
          case 'D': /* left */
            if (_pr.pos>l) {
              _pr.pos--; addch(8);
              if (!prompt_reshow(line)) refresh(); }
            break;
          case 'C': /* right */
            if (_pr.pos<strlen(line)) {
              addch(line[_pr.pos++]);
              prompt_reshow(line); }
            break; } }
    if (c==1) { /* ^a = HOME */
      _pr.point=0, _pr.pos=l;
      prompt_showline(line); }      
    if (c==5) { /* ^e = END */
      _pr.pos=strlen(line);
      prompt_showline(line); }
    if (c==21 || c==25) {
      /* ^u ^y: kill whole line */
      memset(line+l,0,PROMPT_N-l);
      _pr.point=0; _pr.pos=l; 
      prompt_showline(line); }
    if (c==11) {
      /* ^k: delete to eol */
      memset(line+_pr.pos,0,PROMPT_N-_pr.pos);
      prompt_showline(line); }
    if (c==8 || c==127) if (_pr.pos>l) _pr.pos--,c=4;
    if (c==16) ;  /* ^p */
    if (c==4) {
      /* ^d = Delete */
      if (lastc==c && !strcmp(line+l,"exit")) exit(0);
      if (_pr.pos==l && line[l]==0) {
        strcpy(line+l,"exit");
        _pr.pos=l+4; }
      else {
        for (ch=line+_pr.pos; *ch; ch++) *ch=ch[1];
        *ch=ch[l]; }
      prompt_showline(line); }
    
  } while (c!='\n' && c!='\r');

  _pr.back=0;
  _pr.pos=_pr.point=0;
  prompt_showline(line);
  _pr.back=1;
  addch('\n');
  refresh();

  /* the _pr.history */
  if (strcmp(line,_pr.hist[0]) && strcmp(line,_pr.prompt)) {
    for (h=PROMPT_H-1; h>0; h--) strcpy(_pr.hist[h],_pr.hist[h-1]);
    while (line[strlen(line)-1]==' ') line[strlen(line)-1]=0;
    memset(_pr.hist[0],0,PROMPT_N);
    strcpy(_pr.hist[0],line); }

  return line;
}

char *mygetline(char *line) /************************************* mygetline */
{
  char *ch;
  int l=strlen(_pr.prompt);

  prompt_getline(line);
  for (ch=line; *ch; ch++) *ch=ch[l];
  *ch=ch[l];

  return line;
}
