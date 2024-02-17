/* kostnice2 -8 -u < simhelpcz.utf8.c > simhelpcz.c
   for LANG=2 version:
   kostnice2 -8 -a < simhelpcz.utf8.c > simhelpcz0.c
   =============================================
   >>>>>> opravuj soubor simhelpcz.utf8.c <<<<<<
   =============================================
*/
const char *HELP="\
Tento program umo¾òuje studovat skupenské pøemìny látky i dal¹í\n\
molekulové jevy metodami poèítaèových simulací.\n\
\n\
Molekuly se pøitahují, jsou-li tìsnì u sebe.\n\
V na¹em modelu k tomu dojde, kdy¾ se atomy-koleèka témìø dotýkají.\n\
Jsou-li molekuly je¹tì blí¾e u sebe, zaènou se odpuzovat.\n\
V na¹em modelu se odpuzují, pokud se atomy-koleèka pøekrývají.\n\
Jsou-li atomy daleko od sebe, je jejich pøita¾livost ji¾ nepatrná.\n\
\n\
Zkus systém ochlazovat a ohøívat - pozor, chvilku to trvá!\n\
Pevná fáze (krystal) je stabilní okolo teploty 0.05,\n\
kapalina okolo 0.15, plyn okolo 1.\n\
\n\
Teplota (v \"libovolných\" jednotkách, ne ve stupních):\n\
  T = nastavená teplota termostatu\n\
  Tk/Td = teplota systému (kolísá v prùbìhu simulace)\n\
\n\
Graf ukazuje radiální distribuèní funkci, tj. pravdìpodobnost\n\
nalezení dvou èástic v dané vzdálenosti normovanou tak, ¾e pro\n\
èástice bez interakce (ideální plyn) je rovna 1.\n\
\n\
-------------------------------------------------------------------\n\
* Pøehled zkratek a dal¹ích kláves: stiskni h\n\
* Jsou-li pohyby trhané: stiskni a\n\
* Návod ke zvolenému funkènímu tlaèítku získá¹ tak, ¾e na nì\n\
  klikne¹ PRAVÝM my¹ím tlaèítkem.\n\
* Konec: stiskni 2x ESC\n\
-------------------------------------------------------------------\n\
\n\
Teï pokraèuj stiskem jakékoliv klávesy...";

const char *SPEEDHELP="\
Tato tlaèítka ovládají rychlost zobrazování a simulace.\n\
Z klávesnice: s S\n\
\n\
Pøesnìji:\n\
- r=+0 znamená, ¾e se po výpoètu 1 MD kroku nebo 1x hýbnutí èásticí se\n\
  konfigurace zobrazí a vypoètou výsledky (radiální distribuèní funkce\n\
  a tlak)\n\
- r<0 znamená, ¾e se vkládají prodlevy\n\
- r>0, ¾e se nezobrazují v¹echny kroky a výsledky se nepoèítají tak èasto";

const char *THELP="\
Stiskni [ochlaï] nebo otoè my¹í koleèko, chce¹-li sní¾it teplotu.\n\
Stiskni [ohøej] nebo otoè my¹í koleèko, chce¹-li zvý¹it teplotu.\n\
Z klávesnice: t T\n\
\n\
Takto se nastavuje teplota termostatu T. Systému èástic chvilku trvá,\n\
ne¾ se teplota pøiblí¾í k teplotì nastavené.\n\
\n\
Více:\n\
\n\
Teplota je parametrem metody Monte Carlo. Lze si to pøedstavit tak,\n\
¾e pøi jednotlivých pohybech je molekula ve styku s termostatem.\n\
Pøi vy¹¹í teplotì jsou pravdìpodobnìj¹í polohy s vy¹¹í energií, pøi\n\
velmi nízkých teplotách jen ty s nejni¾¹í energií.\n\
\n\
V molekulové dynamice se normálnì zachovává celková (kinetická +\n\
potenciální) energie a teplota je urèena kinetickou energií.\n\
Abychom mohli mít systém o zadané teplotì i pøi simulaci molekulovou\n\
dynamikou, mìníme v ka¾dém kroku integrace o malý kousek v¹echny\n\
rychlosti a tím i teplotu tak, abychom se pøiblí¾ili dané teplotì.\n\
\n\
Stiskem klávesy e se pøepíná na simulaci pøi konstantní energii a zpìt.";

const char *WALLHELP="\
Stiskem tlaèítek se pøepíná pøita¾livá a odpudivá stìna.\n\
Z klávesnice: x X y Y\n\
\n\
Pøita¾livá stìna je zelenì, odpudivá èervenì.\n\
Interakèní parametr stìny lze nastavit v prùbìhu startu.";

const char *VACANCYHELP="\
Ideální krystal s bodovou poruchou - vakancí (chybìjícím atomem\n\
v krystalové møí¾ce). Pøi simulaci porucha cestuje a èasem zmizí na\n\
okraji krystalu.\n\
Z klávesnice: F11\n\
\n\
Nevidí¹-li zøetelnì, kde je porucha, pou¾ij \"uka¾: [sousedy]\".";

const char *NBRHELP="\
Pro ka¾dou molekulu se zobrazí její koordinaèní èíslo, tedy\n\
poèet sousedù do vzdálenosti støedù molekul 1.5 prùmìru molekuly.\n\
Z klávesnice: n\n\
\n\
V dokonalém krystalu je toto èíslo rovno 6 (v reálném trojrozmìrném\n\
systému to je obvykle více, podle krystalové møí¾ky), v kapalinì\n\
kolísá okolo 5 (v trojrozmìrném systému okolo 10-12), v plynu je men¹í\n\
(ve velmi øídkém plynu témìø 0).\n\
\n\
Po chvilièce simulace pokraèuje...\n\
\n\
Pro experty:\n\
Toto èíslo se poèítá z radiální distribuèní funkce, kterou\n\
integrujeme do vzdálenosti 1.5 pøes element rho*2*pi*r*dr.";

const char *METHODHELP="\
K simulaci se pou¾ívaji dvì metody:\n\
\n\
Monte Carlo metoda vzorkuje v¹echny mo¾né polohy molekul pomocí\n\
náhodných pohybù. Náhodné posunutí molekuly je pøijato nebo\n\
odmítnuto podle toho, jaká je energie této nové polohy a jaká je\n\
teplota. Napøíklad pokud by se mìly dvì molekuly pøekrývat, nová\n\
poloha bude odmítnuta.\n\
\n\
Metoda molekulární dynamiky øe¹í pohyb molekul v èase pomocí\n\
numerické integrace Newtonových pohybových rovnic. Molekuly se\n\
pohybují \"jako ve skuteènosti\" - jen mnohokrát zpomalenì.\n\
\n\
Metody se nastavují pøepínaèem [MC/MD] (z klávesnice: m).\n\
Na zaèátku a také po pou¾ití tlaèítek [plyn] atd. je nastaveno\n\
[auto] (z klávesnice: u), co¾ znamená, ¾e vhodná metoda se zvolí\n\
automaticky.\n\
\n\
Velikost posunutí v metodì MC a velikost integraèního kroku se\n\
nastavují automaticky. Toto lze zmìnit stiskem kláves d D.";

const char *GASHELP="\
Stiskem tohoto tlaèítka si vyrobí¹ plyn pøi teplotì T=1.\n\
Z klávesnice: F5\n\
Zkus ochlazovat a pozoruj, k èemu dochází!";

const char *LIQUIDHELP="\
Simulace startuje z malé kapièky. Po ochlazení kapièka\n\
zkrystalizuje, po ohøátí se vypaøí.\n\
Z klávesnice: F6";

const char *CAPILHELP="\
Kapilarita pod mikroskopem.\n\
Z klávesnice: F7\n\
\n\
Dno a postranní stìny jsou pøita¾livé.\n\
Teplota je nastavena na T=0.15, co¾ je v kapalné oblasti.\n\
Dále je pøidána malá tíhová síla.\n\
\n\
Více:\n\
- gravitaci lze ovládat klávesami g G z nebo tlaèítky vedle \"gravitace:\"\n\
- pøita¾livost/odpudivost stìn se pøepíná klávesami x X y Y\n\
  nebo tlaèítky vedle \"stìny:\"\n\
\n\
Tip: padající kapka:\n\
- nastav teplotu na asi 0.13\n\
- horní stìnu udìlej pøita¾livou, ostatní odpudivé\n\
- pomocí kladné gravitace shromá¾di v¹echny molekuly nahoøe\n\
- nastav gravitaci na -0.01: oddìlí se kapièka a spadne";

const char *IDHELP="\
Ideální krystal.\n\
Z klávesnice: F9\n\
Ná¹ dvourozmìrný model krystaluje v ¹estereèné soustavì.\n\
Zkus krystal roztavit - stiskni nìkolikrát tlaèítko [ohøej].";

const char *EDGEHELP="\
Ideální krystal s (dvourozmìrnou) hranovou dislokací - èást jedné\n\
øady atomù chybí. Pøi simulaci porucha cestuje a èasem zmizí na\n\
okraji krystalu.\n\
Z klávesnice: F10\n\
\n\
Nevidí¹-li zøetelnì, kde je porucha, pou¾ij \"uka¾: [sousedy]\".";

const char *INTERHELP="\
Ideální krystal s bodovou poruchou - intersticiálním atomem\n\
(atomem navíc). Pøi simulaci porucha cestuje a èasem zmizí\n\
na okraji krystalu.\n\
Z klávesnice: F12\n\
\n\
Nevidí¹-li zøetelnì, kde je porucha, pou¾ij \"uka¾: [sousedy]\".";

const char *GHELP="\
Tíhová síla (gravitace):\n\
- stiskni [] (nebo klávesu g) pro zvý¹ení gravitace (smìrem dolù)\n\
- stiskni [] (nebo klávesu G) pro sní¾ení gravitace\n\
  (èi zvý¹ení \"antigravitace\" smìrem nahoru)\n\
- Stiskni [0] (nebo klávesu z) pro vynulování gravitace\n\
\n\
Pro plyn (T>1) a zapnutou gravitaci uvidíte, ¾e plyn je\n\
pøi zemi hust¹í ve shodì s barometrickou rovnicí";
