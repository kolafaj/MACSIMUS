char *Mendeleyev(double m)
/* atom symbol from the mass */
{
  static struct Mendeleyev_s {
    char *At;
    int Z;
    double M; } M[]={
    {"H",        1,        1.00794},    
    {"D",        1,        2.0141017},   
    {"He",       2,        4.002602},    
    {"Li",       3,        6.941},      
    {"Be",       4,        9.012182},    
    {"B",        5,       10.811},       
    {"N",        7,       14.00674},     
    {"O",        8,       15.9994},      
    {"F",        9,       18.9984032},   
    {"Ne",      10,       20.1797},     
    {"Na",      11,       22.989768},   
    {"Mg",      12,      24.3050},      
    {"Al",      13,      26.981539},    
    {"Si",      14,      28.0855},      
    {"P",       15,      30.973762},    
    {"S",       16,       32.066},      
    {"Cl",      17,      35.4527},      
    {"Ar",      18,      39.948},       
    {"K",       19,      39.0983},      
    {"Ca",      20,       40.078},      
    {"Sc",      21,       44.955910},   
    {"Ti",      22,       47.88},       
    {"V",       23,       50.9415},     
    {"Cr",      24,      51.9961},      
    {"Mn",      25,       54.93805},    
    {"Fe",      26,       55.847},
    {"Co",      27,      58.93320},     
    {"Ni",      28,       58.6934},     
    {"Cu",      29,       63.546},      
    {"Zn",      30,       65.39},       
    {"Ga",      31,      69.723},       
    {"Ge",      32,      76.61},        
    {"As",      33,      74.92159},     
    {"Se",      34,       78.96},       
    {"Br",      35,      79.904},       
    {"Kr",      36,      83.80},        
    {"Rb",      37,       85.4678},     
    {"Sr",      38,       87.62},       
    {"Y",       39,      88.90585},     
    {"Zr",      40,       91.224},      
    {"Nb",      41,       92.90638},    
    {"Mo",      42,       95.94},       
    {"Tc",      43,       97.9072},     
    {"Ru",      44,      101.07},       
    {"Rh",      45,      102.90550},    
    {"Pd",      46,      106.42},       
    {"Ag",      47,      107.8682},     
    {"Cd",      48,     112.411},       
    {"In",      49,     114.818},       
    {"Sn",      50,     118.710},       
    {"Sb",      51,     121.757},       
    {"Te",      52,      127.60},       
    {"I",       53,     126.90447},     
    {"Xe",      54,      131.29},       
    {"Cs",      55,     132.90543},     
    {"Ba",      56,     137.327},       
    {"La",      57,      138.9055},     
    {"Ce",      58,     140.115},       
    {"Pr",      59,      140.90765},    
    {"C",       6,        12.011},      
    {"Nd",      60,      144.24},       
    {"Pm",      61,      144.9127},     
    {"Sm",      62,      150.36},       
    {"Eu",      63,     151.965},       
    {"Gd",      64,     157.25},        
    {"Tb",      65,      158.92534},    
    {"Dy",      66,     162.50},        
    {"Ho",      67,     164.93032},     
    {"Er",      68,     167.26},        
    {"Tm",      69,      168.93421},    
    {"Yb",      70,      173.04},       
    {"Lu",      71,      174.967},      
    {"Hf",      72,     178.49},        
    {"Ta",      73,      180.9479},     
    {"W",       74,      183.84},       
    {"Re",      75,      186.207},      
    {"Os",      76,      190.23},       
    {"Ir",      77,     192.22},        
    {"Pt",      78,      195.08},       
    {"Au",      79,      196.96654},    
    {"Hg",      80,      200.59},       
    {"Tl",      81,      204.3833},     
    {"Pb",      82,      207.2},        
    {"Bi",      83,     208.98037},     
    {"Po",      84,      208.9824},     
    {"At",      85,     209.9871},      
//{"Fr",        87, 0 },        
    {"Ra",      88,      226.0254},     
    {"Ac",      89,     227.0278},      
    {"Th",      90,      232.0381},     
    {"Pa",      91,      231.03588},    
    {"U",       92,      238.0289},     
    {"Np",      93,      237.0482},     
    {"Pu",      94,      244.0642},     
    {"Am",      95,     243.0614},      
    {"Cm",      96,     247.0703},      
    {"Bk",      97,     247.0703},      
    {"Cf",      98,     251.0796},      
    {"Es",      99,     252.083},       
    {"Fm",      100 ,   257.0951},      
    {"Md",      101,     258.10},       
    {"No",      102,     259.1009},     
    {"Lr",      103,     262.11},       
    {"Unq",     104,     261.11},       
    {"Unp",     105,     262.114},      
    {"Unh",     106,     263.118},      
    {"Uns",     107,     262.12},       
    {"CH1",0,13.019},
    {"CH2",0,14.027},
    {"CH3",0,15.035},
    {NULL,      0,     0}};
      //{"Rn",  86, 0 },        
  int i,ii;
  double d=999;

  for (i=0; M[i].At; i++)
    if (fabs(M[i].M-m)<d) ii=i,d=fabs(M[i].M-m);

  return M[ii].At;
}
