#ifdef DEBUG
prt("DEBUG NB: E=En.pot+En.el+En.self, En.vir is without elst energy");
    prt("DEBUG V          =%16.6f",box.V);
    prt("DEBUG Punit/3/V  =%16.6f",Punit/3/box.V);
    //    prt("DEBUG r[0]       =%10.6f %10.6f %10.6f",VARG(rof(cfg,cfg[0]->rp)[0]));
    prt("DEBUG dE/dV      =%16.6f",(Ep-Em)/(2*constrd.dV));
    prt("DEBUG dEn.el/dV  =%16.6f",(Enp.el-En.el)/(2*constrd.dV));
    prt("DEBUG dEn.pot/dV =%16.6f",(Enp.pot-En.pot)/(2*constrd.dV));
    //    prt("DEBUG Ep        =%16.6f",Ep);
    //    prt("DEBUG Em        =%16.6f",Em);
    prt("DEBUG av.E       =%16.6f",(Em+Ep)/2);
    //    prt("DEBUG Enp.pot   =%16.6f",Enp.pot);
    //    prt("DEBUG En.pot    =%16.6f",En.pot);
    prt("DEBUG av.En.pot  =%16.6f",(En.pot+Enp.pot)/2);
    //    prt("DEBUG Enp.el    =%16.6f",Enp.el);
    //    prt("DEBUG En.el     =%16.6f",En.el);
    prt("DEBUG av.En.el   =%16.6f",(En.el+Enp.el)/2);
    prt("DEBUG av.En.vir  =%16.6f",(En.vir+Enp.vir)/2);
    //    prt("DEBUG av.En.sum  =%16.6f",(En.el+En.pot+En.self+Enp.el+Enp.pot+Enp.self)/2);
#  ifdef POLAR
    prt("DEBUG dEn.self/dV=%16.6f",(Enp.self-En.self)/(2*constrd.dV));
 /*.....   prt("DEBUG both*3:=%16.6f",-(En.el+En.self-Enp.self-Enp.el)/(2*constrd.dV)*3);*/
    prt("DEBUG av.En.self =%16.6f",(En.self+Enp.self)/2);
#  endif /*# POLAR */
#endif /*# DEBUG */
