
#include <stdio.h>

/*
 * get option letter from argument vector
 */
int	opt_err = 1,		/* useless, never set or used */
	opt_ind = 1,		/* index into parent argv vector */
	opt_opt;			/* character checked for validity */
char	*opt_arg;		/* argument associated with option */

#define BADCH	(int)'?'
#define EMSG	""
#define tell(s)	fputs(*nargv,stderr);fputs(s,stderr); \
		fputc(opt_opt,stderr);fputc('\n',stderr);return(BADCH);

get_opt(nargc,nargv,ostr)
int	nargc;
char	**nargv,
	*ostr;
{
	static char	*place = EMSG;	/* option letter processing */
	register char	*oli;		/* option letter list index */
	char	*index();

	if(!*place) {			/* update scanning pointer */
		if(opt_ind >= nargc || *(place = nargv[opt_ind]) != '-' || !*++place) return(EOF);
		if (*place == '-') {	/* found "--" */
			++opt_ind;
			return(EOF);
		}
	}				/* option letter okay? */
	if ((opt_opt = (int)*place++) == (int)':' || !(oli = index(ostr,opt_opt))) {
		if(!*place) ++opt_ind;
		tell(": illegal option -- ");
	}
	if (*++oli != ':') {		/* don't need argument */
		opt_arg = NULL;
		if (!*place) ++opt_ind;
	}
	else {				/* need an argument */
		if (*place) opt_arg = place;	/* no white space */
		else if (nargc <= ++opt_ind) {	/* no arg */
			place = EMSG;
			tell(": option requires an argument -- ");
		}
	 	else opt_arg = nargv[opt_ind];	/* white space */
		place = EMSG;
		++opt_ind;
	}
	return(opt_opt);			/* dump back option letter */
}
