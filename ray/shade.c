/***********************************************************************
 * $Author: markv $
 * $Revision: 1.2 $
 * $Date: 88/09/12 12:58:01 $
 * $Log:        shade.c,v $
 * Revision 1.2  88/09/12  12:58:01  markv
 * Added specular reflections and shadow caching.  When a object
 * is found between a light source and the current point we are trying
 * to shade, that object is cached (indexed by recursion level) in the
 * lightsource.  Next time we test a shadow against the light source,
 * we test this object first.  If it is between us and the light source,
 * we can correctly shadow the object without calling Intersect().
 * 
 * Note: specular highlights call the unix pow() function, which seems
 * to be REALLY expensive.  Optimizations could be made here.
 * 
 * Revision 1.1  88/09/11  11:00:44  markv
 * Initial revision
 *  
 ***********************************************************************/

#include <stdio.h>
#include <math.h>
#include "defs.h"
#include "extern.h"

/***********************************************************************
 * Shade(level, weight, P, N, I, hit, col)
 *
 * Wow! Too many parameters!
 *
 * Shade is the driving force of the raytracer.  It calculates the actual
 * luminance at point P, given N, the normal, and I the incident vector.
 * We also pass in the hit, because the normal generating code might need it.
 * The color is returned in col.
 ***********************************************************************/

Shade(level, weight, P, N, I, hit, col) 
 int level ;
 Flt weight ;
 Vec P, N, I ;
 Isect * hit;
 Color col ;
{
        Ray     tray ;
        Color   tcol ;
        Vec     L, /*H,*/ R ;
        Flt     t ;
        Flt     diff ;
        Flt     spec ;
#ifdef SHADOW_CACHING
        Object  * cached ;
#endif /* SHADOW_CACHING */
        Isect nhit ;
        Surface * surf ;
        int l ;

/*.....        col[0] = col[1] = col[2] = 0.0 ;*/

        surf = hit -> isect_prim -> o_surf ;

        { /* JK: ambient light */
          double f=fabs(I[0]*N[0]+I[1]*N[1]+I[2]*N[2])
            /sqrt((I[0]*I[0]+I[1]*I[1]+I[2]*I[2])
                  *(N[0]*N[0]+N[1]*N[1]+N[2]*N[2]))*normaldirlight
	    +isotropiclight;
          for (l=0; l<3; l++) col[l]=surf->surf_color[l]*f;
        }

        for (l = 0; l < nLights; l++) {
                /* JK: */
                if (diffuselight!=0) {
                  Vec r;
                  int i;
                  for (i=0; i<3; i++) {
                    r[i]=(rndcos()-rndcos()+rndcos()-rndcos()+rndcos()-rndcos())*diffuselight;
                    VecAdd(Lights[l].light_pos0,r,Lights[l].light_pos); }
                }
                VecSub(Lights[l].light_pos, P, L);
                if (VecDot(N,L) >= 0.0) {
                        t = VecNormalize(L);
                        VecCopy(P, tray.P);
                        VecCopy(L, tray.D);
                        nShadows ++ ;
#ifdef SHADOW_CACHING
                        cached = Lights[l].light_obj_cache[level] ;
                        if (cached 
                                && (cached -> o_procs -> intersect) 
                                        (cached, &tray, &nhit) 
                                && nhit.isect_t < t) {
                                /* 
                                 * we are in shadow, continue...
                                 */
                                nShadowCacheHits ++ ;
                                continue ;
                        }
#endif /* SHADOW_CACHING */
                        if (Shadow(&tray, &nhit, t)) {
                                diff = VecDot(N,L) * surf -> surf_kd 
                                        * Lights[l].light_brightness ;
#ifdef SHADOW_CACHING
                                Lights[l].light_obj_cache[level] = NULL ;
#endif /* SHADOW_CACHING */

#if 1 /* normal */
                                VecAddS(diff, surf -> surf_color, col, col) ;
#else /* very special patch */
if (surf -> surf_color[1]>1.009) {
  Color ccc;
  int i;
  for (i=0; i<3; i++) ccc[i]=surf -> surf_color[i]*0;
  ccc[0]=0.5+(cos(9.*tray.P[0])*cos(11.*tray.P[2]))*0.7;
  VecAddS(diff, ccc, col, col) ; }
else
  VecAddS(diff, surf -> surf_color, col, col) ;
#endif

                                SpecularDirection(I, N, R) ;
                                VecNormalize(R) ;
                                if (surf -> surf_shine > rayeps) {
                                        spec = VecDot(R,L) ;
                                        if (spec > rayeps) {
                                                spec = pow(spec, surf -> surf_shine ) * Lights[l].light_brightness ;
                                                col[0] += spec ; 
                                                col[1] += spec ; 
                                                col[2] += spec ;
                                        }
                                }
                        } else {
#ifdef SHADOW_CACHING
                                Lights[l].light_obj_cache[level] = 
                                        nhit.isect_prim ;
#endif /* SHADOW_CACHING */
                        }
                                
                }
        }

	/* fog - not very good, suitable for background only */
	if (fogthick) {
	  double f=P[2]-fog;

	  if (f>0) f=1;
	  else f=exp(f/fogthick);

	  for (l=0; l<3; l++) col[l]=col[l]*f+(1-f)*BackgroundColor[l];
	}

        VecCopy(P, tray.P);

        if(surf -> surf_ks * weight > minweight) {
                nReflected ++ ;
                SpecularDirection(I, N, tray.D);
                VecNormalize(tray.D);
                Trace(level + 1, surf -> surf_ks * weight, &tray, tcol);
                VecAddS(surf -> surf_ks, tcol, col, col);
                
        }

#if 0 /* not very successful attempt to handle diffuse light */
        if (level==0) { /* JK: */
                int i; Vec r; Ray ttray=tray;
                SpecularDirection(I, N, tray.D);
                VecNormalize(tray.D);
#define NX 20
                for (i=0; i<NX; i++) {
                  r[0]=rndcos()-rndcos();
                  r[1]=rndcos()-rndcos();
                  r[2]=rndcos()-rndcos();
                  VecAddS(0.3,tray.D,r,ttray.D);
                  VecNormalize(ttray.D);
                  Trace(level + 1, weight, &ttray, tcol);
                  VecAddS(0.5 /*?*/ /NX, tcol, col, col); } }
#endif

        if (surf -> surf_kt * weight > minweight) { 
                nRefracted ++ ;
                if (hit -> isect_enter) 
                        TransmissionDirection(NULL, surf, I, N, tray.D) ;
                else    
                        TransmissionDirection(surf, NULL, I, N, tray.D) ;
                Trace(level + 1, surf -> surf_kt * weight, &tray, tcol) ;
                VecAddS(surf -> surf_kt, tcol, col, col) ;
        }
return 0;
}

/***********************************************************************
 * SpecularDirection(I, N, R)
 * 
 * Given an incident vector I, and the normal N, calculate the 
 * direction of the reflected ray R.
 ***********************************************************************/

int
SpecularDirection(I, N, R)
 Vec I, N, R;
{
        VecComb(1.0/fabs(VecDot(I,N)), I, 2.0, N, R);
        VecNormalize(R);
        return 0;
}

/***********************************************************************
 * TransmissionDirection(m1, m2, I, N, T)
 *
 * calculates the direction of the transmitted ray
 ***********************************************************************/

TransmissionDirection(m1, m2, I, N, T)
 Surface *m1, *m2;
 Vec I, N, T ;
{
        Flt n1, n2, eta, c1, cs2 ;
        n1 = m1 ? m1 -> surf_ior : 1.0 ;
        n2 = m2 ? m2 -> surf_ior : 1.0 ;
        eta = n1/n2 ;

        c1 = -VecDot(I,N);
        cs2 = 1.0 - eta * eta*(1.0 - c1*c1);
        if (cs2 < 0.0)
                return 0;
        VecComb(eta, I, eta*c1-sqrt(cs2), N, T);
        return(1);
}
