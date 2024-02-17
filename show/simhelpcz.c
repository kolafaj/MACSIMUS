/* kostnice2 -8 -u < simhelpcz.utf8.c > simhelpcz.c
   for LANG=2 version:
   kostnice2 -8 -a < simhelpcz.utf8.c > simhelpcz0.c
   =============================================
   >>>>>> opravuj soubor simhelpcz.utf8.c <<<<<<
   =============================================
*/
const char *HELP="\
Tento program umo��uje studovat skupensk� p�em�ny l�tky i dal��\n\
molekulov� jevy metodami po��ta�ov�ch simulac�.\n\
\n\
Molekuly se p�itahuj�, jsou-li t�sn� u sebe.\n\
V na�em modelu k tomu dojde, kdy� se atomy-kole�ka t�m�� dot�kaj�.\n\
Jsou-li molekuly je�t� bl�e u sebe, za�nou se odpuzovat.\n\
V na�em modelu se odpuzuj�, pokud se atomy-kole�ka p�ekr�vaj�.\n\
Jsou-li atomy daleko od sebe, je jejich p�ita�livost ji� nepatrn�.\n\
\n\
Zkus syst�m ochlazovat a oh��vat - pozor, chvilku to trv�!\n\
Pevn� f�ze (krystal) je stabiln� okolo teploty 0.05,\n\
kapalina okolo 0.15, plyn okolo 1.\n\
\n\
Teplota (v \"libovoln�ch\" jednotk�ch, ne ve stupn�ch):\n\
  T = nastaven� teplota termostatu\n\
  Tk/Td = teplota syst�mu (kol�s� v pr�b�hu simulace)\n\
\n\
Graf ukazuje radi�ln� distribu�n� funkci, tj. pravd�podobnost\n\
nalezen� dvou ��stic v dan� vzd�lenosti normovanou tak, �e pro\n\
��stice bez interakce (ide�ln� plyn) je rovna 1.\n\
\n\
-------------------------------------------------------------------\n\
* P�ehled zkratek a dal��ch kl�ves: stiskni h\n\
* Jsou-li pohyby trhan�: stiskni a\n\
* N�vod ke zvolen�mu funk�n�mu tla��tku z�sk� tak, �e na n�\n\
  klikne� PRAV�M my��m tla��tkem.\n\
* Konec: stiskni 2x ESC\n\
-------------------------------------------------------------------\n\
\n\
Te� pokra�uj stiskem jak�koliv kl�vesy...";

const char *SPEEDHELP="\
Tato tla��tka ovl�daj� rychlost zobrazov�n� a simulace.\n\
Z kl�vesnice: s S\n\
\n\
P�esn�ji:\n\
- r=+0 znamen�, �e se po v�po�tu 1 MD kroku nebo 1x h�bnut� ��stic� se\n\
  konfigurace zobraz� a vypo�tou v�sledky (radi�ln� distribu�n� funkce\n\
  a tlak)\n\
- r<0 znamen�, �e se vkl�daj� prodlevy\n\
- r>0, �e se nezobrazuj� v�echny kroky a v�sledky se nepo��taj� tak �asto";

const char *THELP="\
Stiskni [ochla�] nebo oto� my�� kole�ko, chce�-li sn�it teplotu.\n\
Stiskni [oh�ej] nebo oto� my�� kole�ko, chce�-li zv��it teplotu.\n\
Z kl�vesnice: t T\n\
\n\
Takto se nastavuje teplota termostatu T. Syst�mu ��stic chvilku trv�,\n\
ne� se teplota p�ibl�� k teplot� nastaven�.\n\
\n\
V�ce:\n\
\n\
Teplota je parametrem metody Monte Carlo. Lze si to p�edstavit tak,\n\
�e p�i jednotliv�ch pohybech je molekula ve styku s termostatem.\n\
P�i vy��� teplot� jsou pravd�podobn�j�� polohy s vy��� energi�, p�i\n\
velmi n�zk�ch teplot�ch jen ty s nejni��� energi�.\n\
\n\
V molekulov� dynamice se norm�ln� zachov�v� celkov� (kinetick� +\n\
potenci�ln�) energie a teplota je ur�ena kinetickou energi�.\n\
Abychom mohli m�t syst�m o zadan� teplot� i p�i simulaci molekulovou\n\
dynamikou, m�n�me v ka�d�m kroku integrace o mal� kousek v�echny\n\
rychlosti a t�m i teplotu tak, abychom se p�ibl�ili dan� teplot�.\n\
\n\
Stiskem kl�vesy e se p�ep�n� na simulaci p�i konstantn� energii a zp�t.";

const char *WALLHELP="\
Stiskem tla��tek se p�ep�n� p�ita�liv� a odpudiv� st�na.\n\
Z kl�vesnice: x X y Y\n\
\n\
P�ita�liv� st�na je zelen�, odpudiv� �erven�.\n\
Interak�n� parametr st�ny lze nastavit v pr�b�hu startu.";

const char *VACANCYHELP="\
Ide�ln� krystal s bodovou poruchou - vakanc� (chyb�j�c�m atomem\n\
v krystalov� m��ce). P�i simulaci porucha cestuje a �asem zmiz� na\n\
okraji krystalu.\n\
Z kl�vesnice: F11\n\
\n\
Nevid�-li z�eteln�, kde je porucha, pou�ij \"uka�: [sousedy]\".";

const char *NBRHELP="\
Pro ka�dou molekulu se zobraz� jej� koordina�n� ��slo, tedy\n\
po�et soused� do vzd�lenosti st�ed� molekul 1.5 pr�m�ru molekuly.\n\
Z kl�vesnice: n\n\
\n\
V dokonal�m krystalu je toto ��slo rovno 6 (v re�ln�m trojrozm�rn�m\n\
syst�mu to je obvykle v�ce, podle krystalov� m��ky), v kapalin�\n\
kol�s� okolo 5 (v trojrozm�rn�m syst�mu okolo 10-12), v plynu je men��\n\
(ve velmi ��dk�m plynu t�m�� 0).\n\
\n\
Po chvili�ce simulace pokra�uje...\n\
\n\
Pro experty:\n\
Toto ��slo se po��t� z radi�ln� distribu�n� funkce, kterou\n\
integrujeme do vzd�lenosti 1.5 p�es element rho*2*pi*r*dr.";

const char *METHODHELP="\
K simulaci se pou��vaji dv� metody:\n\
\n\
Monte Carlo metoda vzorkuje v�echny mo�n� polohy molekul pomoc�\n\
n�hodn�ch pohyb�. N�hodn� posunut� molekuly je p�ijato nebo\n\
odm�tnuto podle toho, jak� je energie t�to nov� polohy a jak� je\n\
teplota. Nap��klad pokud by se m�ly dv� molekuly p�ekr�vat, nov�\n\
poloha bude odm�tnuta.\n\
\n\
Metoda molekul�rn� dynamiky �e�� pohyb molekul v �ase pomoc�\n\
numerick� integrace Newtonov�ch pohybov�ch rovnic. Molekuly se\n\
pohybuj� \"jako ve skute�nosti\" - jen mnohokr�t zpomalen�.\n\
\n\
Metody se nastavuj� p�ep�na�em [MC/MD] (z kl�vesnice: m).\n\
Na za��tku a tak� po pou�it� tla��tek [plyn] atd. je nastaveno\n\
[auto] (z kl�vesnice: u), co� znamen�, �e vhodn� metoda se zvol�\n\
automaticky.\n\
\n\
Velikost posunut� v metod� MC a velikost integra�n�ho kroku se\n\
nastavuj� automaticky. Toto lze zm�nit stiskem kl�ves d D.";

const char *GASHELP="\
Stiskem tohoto tla��tka si vyrob� plyn p�i teplot� T=1.\n\
Z kl�vesnice: F5\n\
Zkus ochlazovat a pozoruj, k �emu doch�z�!";

const char *LIQUIDHELP="\
Simulace startuje z mal� kapi�ky. Po ochlazen� kapi�ka\n\
zkrystalizuje, po oh��t� se vypa��.\n\
Z kl�vesnice: F6";

const char *CAPILHELP="\
Kapilarita pod mikroskopem.\n\
Z kl�vesnice: F7\n\
\n\
Dno a postrann� st�ny jsou p�ita�liv�.\n\
Teplota je nastavena na T=0.15, co� je v kapaln� oblasti.\n\
D�le je p�id�na mal� t�hov� s�la.\n\
\n\
V�ce:\n\
- gravitaci lze ovl�dat kl�vesami g G z nebo tla��tky vedle \"gravitace:\"\n\
- p�ita�livost/odpudivost st�n se p�ep�n� kl�vesami x X y Y\n\
  nebo tla��tky vedle \"st�ny:\"\n\
\n\
Tip: padaj�c� kapka:\n\
- nastav teplotu na asi 0.13\n\
- horn� st�nu ud�lej p�ita�livou, ostatn� odpudiv�\n\
- pomoc� kladn� gravitace shrom�di v�echny molekuly naho�e\n\
- nastav gravitaci na -0.01: odd�l� se kapi�ka a spadne";

const char *IDHELP="\
Ide�ln� krystal.\n\
Z kl�vesnice: F9\n\
N� dvourozm�rn� model krystaluje v �estere�n� soustav�.\n\
Zkus krystal roztavit - stiskni n�kolikr�t tla��tko [oh�ej].";

const char *EDGEHELP="\
Ide�ln� krystal s (dvourozm�rnou) hranovou dislokac� - ��st jedn�\n\
�ady atom� chyb�. P�i simulaci porucha cestuje a �asem zmiz� na\n\
okraji krystalu.\n\
Z kl�vesnice: F10\n\
\n\
Nevid�-li z�eteln�, kde je porucha, pou�ij \"uka�: [sousedy]\".";

const char *INTERHELP="\
Ide�ln� krystal s bodovou poruchou - interstici�ln�m atomem\n\
(atomem nav�c). P�i simulaci porucha cestuje a �asem zmiz�\n\
na okraji krystalu.\n\
Z kl�vesnice: F12\n\
\n\
Nevid�-li z�eteln�, kde je porucha, pou�ij \"uka�: [sousedy]\".";

const char *GHELP="\
T�hov� s�la (gravitace):\n\
- stiskni [] (nebo kl�vesu g) pro zv��en� gravitace (sm�rem dol�)\n\
- stiskni [] (nebo kl�vesu G) pro sn�en� gravitace\n\
  (�i zv��en� \"antigravitace\" sm�rem nahoru)\n\
- Stiskni [0] (nebo kl�vesu z) pro vynulov�n� gravitace\n\
\n\
Pro plyn (T>1) a zapnutou gravitaci uvid�te, �e plyn je\n\
p�i zemi hust�� ve shod� s barometrickou rovnic�";
