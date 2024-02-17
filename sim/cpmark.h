/* legacy: IEEE endian-reversible, bytes=253 83 83 253, string="ıSSı" */
/* legacy: 1st line and ASCII date/time, NCP>=7 */
#define CPmark ((float)-1.755645e+37)

/* mark time_t (long unsigned) time, bytes=253 84 84 253, string="ıTTı" */
#define CPmarkT ((float)-1.76398512e+37)

/* mark double time t (float for NCP=2), bytes=253 85 85 253, string="ıUUı" */
#define CPmarkU ((float)-1.77232525e+37)

/* mark double cycle noint*h (float for NCP=2), bytes=253 86 86 253, string="ıVVı" */
#define CPmarkV ((float)-1.78066538e+37)

/* mark NCP>=5: both t and h*noint (double) in one line bytes=253 87 87 253, string="ıWWı" */
#define CPmarkW ((float)-1.7890055e+37)

#if 0
/* not used, bytes=253 88 88 253, string="ıXXı" */
#define CPmarkX ((float)-1.79734563e+37)
/* not used, bytes=253 89 89 253, string="ıYYı" */
#define CPmarkY ((float)-1.79734563e+37)
/* not used, bytes=253 90 90 253, string="ıZZı" */
#define CPmarkZ ((float)-1.81402588e+37)
#endif

/* #define VeryOldCPmark -7.654e37 */
