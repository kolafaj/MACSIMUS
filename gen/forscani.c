  /* #included by forscanf.c */
  va_list args;
  char *fmt=(char*)format;
  int c,j,i=0,len,eol=0,l;

  char field[82],*f;

  va_start(args,format);

 again:

  switch (*fmt) {

    case 0: goto done;

    case '\n':
      if (eol) { eol=0; break; /* already eol */ }
      do MYGETC0 while (c!='\n');
      eol=0;
      break;

    case ' ':   
      MYGETC
      break;

    case '%':
      len=atoi(++fmt);
      if (len<0 || len>80) { i=-3; goto done; }
      while (*fmt>='0' && *fmt<='9') fmt++;

      for (j=0; j<len; j++) {
        MYGETC
        field[j]=c; }
      field[j]=0;

      l=0;
      if (*fmt=='l') { l=1; fmt++; } /* arg is long */

      switch (*fmt) {

        case 'd': case 'i':
          if (l) *va_arg(args,long*)=atol(field);
          else *va_arg(args,int*)=atoi(field);
          break;

        case 'f':
          if (l) *va_arg(args,double*)=atof(field);
          else *va_arg(args,float*)=atof(field);
          break;
          
        case 'a':
          memcpy(va_arg(args,char*),field,len);
          break;

        case 's':
          f=field+len-1;
          while (*f==' ' && f>=field) *f--=0;
          f=field;
          while (*f && *f==' ') f++;
          strcpy(va_arg(args,char*),f);
          break;

        default:
          i=-4; goto done; }

      i++;
      break;

    default:
      i=-5; goto done;
	
  }

  fmt++;
  goto again;

 done:
  va_end(args);

  return i;
