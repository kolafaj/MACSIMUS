/*** `?'-statements */
if (*_Id.b=='?') {
  _Id.b++;

  if (*_Id.b<' ') { /*** ? statement single on a line = help ***/
    prt("=== getdata/calc V.%d ===",CALC);
    prts(
"EXPR     evaluate expression in usual mathematical syntax\n"
"         can use known not-subscripted identifiers and pi\n"
"         DD::HH:MM:SS=time in hours, 0x##=hex. number, 0##=octal\n"
"         ^=**=power, |=or, &=and, $=xor, %=modulo\n"
"ID=EXPR  assign EXPR to ID\n"
"?ID=EXPR assign with echo\n"
"ID+=EXPR add to; similarly -= *= /= %= ^=\n"
";        end data (program will continue)\n"
"ID       print ID value\n"
"#        last value; example: a=exp(#)\n"
"\?\?       list all ID\'s\n\
?=       toggle auto echo\n\
?%FMT  	 set format, e.g. ?%9.7f\n\
!        comment to end of line");

#ifdef SCR
    if (_S_buf) prts("$?       help of scroll");
#endif

    prts("  example:\nid1=2.3*cos(PI/3) id2 += 1 id3=id1^3;");
    prtecho();
    goto noscan; }

  else if (*_Id.b=='=') { /*** ?= statement */
    _Id.b++;
    _Id.echo=(_Id.echo+1)%3;
    prtecho();
    goto noscan; }

  else if (*_Id.b=='%') { /*** ?% statement */
    for (x=_Id.b; *x>' ' && x-_Id.b<13; x++) if (strchr("eEgGf",*x)) {
      x=_Id.fmt+3;
      do *x++=*_Id.b++; while (*_Id.b>' ');
      *x=0;
      if (x-_Id.fmt>=16) Error("format too long");
      goto noscan; }
    Error("bad format");
    *_Id.b=0;
    goto noscan; }

  else if (*_Id.b=='?') { /*** ?? statement */
    _Id.b++;
    _Id.prt=2;     /* mark the list-of-ID function */
    _Id.used=1;    /* suppress check `identifier found' */
    _Id.col=0;     /* (=20000 == pretend eoln - will print \n) */
    goto ret; }

  /* ?<anything else> is accepted for compatibility reasons */
  while (*_Id.b==' ') _Id.b++;
  _Id.id=_Id.b;
  _Id.prt=1; /* print result */

} /* end of ?-statements */
