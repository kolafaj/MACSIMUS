/* the same as strtod,atof but accepting also FORTRAN numbers as 1.2D+3 */
double fstrtod(const char *nptr, char **endptr);
double fatof(const char *nptr);
