#if PRESSURETENSOR&PT_VIR
        if (measure) {
          /* WARNING: with ~h^2 error */
          PTvirc[0]+=(sq/bondq)*(rij[0]*rij[0]); 
          PTvirc[1]+=(sq/bondq)*(rij[1]*rij[1]); 
          PTvirc[2]+=(sq/bondq)*(rij[2]*rij[2]); 
#if PRESSURETENSOR&PT_OFF
          PTvirc[3]+=(sq/bondq)*(rij[1]*rij[2]);
          PTvirc[4]+=(sq/bondq)*(rij[2]*rij[0]);
          PTvirc[5]+=(sq/bondq)*(rij[0]*rij[1]);
#endif
          En.virc += sq; }
#else
        En.virc += sq; /* not needed if !measure, but faster to sum anyway */
#endif        
