/* 
  8x8 font from Turbo Pascal 5.5, in columns, bit0=top (low y) 
  8x9 version: gjpqy with one pixel below line
  V modified
  see outtext0.c or use bgifont.pas
  437 codepage coding
  bytes or shorts are COLUMNS, lsb on top
*/
static unsigned short font_fbgi[256*8]={
0,0,0,0,0,0,0,0,
60,66,149,177,177,149,66,60,
60,126,235,207,207,235,126,60,
14,31,62,124,62,31,14,0,
8,28,62,127,62,28,8,0,
28,28,74,127,127,74,28,28,
8,28,94,127,127,94,28,8,
0,0,60,60,60,60,0,0,
255,255,195,195,195,195,255,255,
0,126,66,66,66,66,126,0,
255,129,189,189,189,189,129,255,
112,248,136,136,249,127,7,15,
78,95,241,241,95,78,0,0,
192,224,254,127,5,5,7,7,
192,254,127,5,5,101,127,63,
219,219,60,231,231,60,219,219,
127,127,62,62,28,8,8,0,
8,8,28,62,62,127,127,0,
36,102,255,255,102,36,0,0,
0,95,95,0,95,95,0,0,
14,31,17,127,127,1,127,127,
0,26,191,165,165,253,88,0,
112,112,112,112,112,112,112,0,
0,148,182,255,255,182,148,0,
4,6,127,127,6,4,0,0,
16,48,127,127,48,16,0,0,
8,8,8,8,62,28,8,0,
8,28,62,8,8,8,8,0,
28,28,16,16,16,16,16,0,
8,28,42,8,8,42,28,8,
48,56,60,62,60,56,48,0,
6,14,30,62,30,14,6,0,
0,0,0,0,0,0,0,0,
0,0,0,95,95,0,0,0,
7,7,0,0,7,7,0,0,
84,126,63,85,126,63,21,0,
36,46,106,107,43,58,18,0,
67,99,48,24,12,102,99,0,
48,122,79,93,55,114,80,0,
0,4,7,3,0,0,0,0,
0,0,28,62,99,65,0,0,
0,0,65,99,62,28,0,0,
8,42,62,28,62,42,8,0,
0,8,8,62,62,8,8,0,
0,0,128,224,96,0,0,0,
0,8,8,8,8,8,8,0,
0,0,0,96,96,0,0,0,
64,96,48,24,12,6,3,0,
62,127,65,73,65,127,62,0, /* 0 : central . added 11/2090*/
0,64,66,127,127,64,64,0,
98,115,81,89,73,111,102,0,
34,99,65,73,73,127,54,0,
24,28,22,19,127,127,16,0,
39,103,69,69,69,125,57,0,
60,126,75,73,73,121,48,0,
3,3,113,121,13,7,3,0,
54,127,73,73,73,127,54,0,
6,79,73,73,105,63,30,0,
0,0,0,102,102,0,0,0,
0,0,128,230,102,0,0,0,
0,8,28,54,99,65,0,0,
0,36,36,36,36,36,36,0,
0,0,65,99,54,28,8,0,
0,2,3,81,89,15,6,0,
62,127,65,93,93,95,30,0,
124,126,19,17,19,126,124,0,
127,127,73,73,73,127,54,0,
62,127,65,65,65,103,38,0,
127,127,65,65,99,62,28,0,
127,127,73,73,73,73,65,0,
127,127,9,9,9,9,1,0,
62,127,65,65,73,123,122,0,
127,127,8,8,8,127,127,0,
0,65,127,127,65,0,0,0,
48,112,64,65,65,127,63,0,
127,127,8,28,54,99,65,0,
127,127,64,64,64,64,64,0,
127,127,6,12,6,127,127,0,
127,127,6,12,24,127,127,0,
62,127,65,65,65,127,62,0,
127,127,9,9,9,15,6,0,
62,127,65,65,65,255,190,0,
127,127,9,9,9,127,118,0,
34,103,77,89,115,34,0,0,
1,1,127,127,1,1,0,0,
63,127,64,64,64,127,63,0,
 /*.....{31,63,96,64,96,63,31,0}, V */
15,31,48,96,48,31,15,0,
127,127,48,24,48,127,127,0,
99,119,28,8,28,119,99,0,
3,7,12,120,120,12,7,3,
97,113,89,77,71,67,65,0,
0,0,127,127,65,65,0,0,
1,3,6,12,24,48,96,64,
0,0,65,65,127,127,0,0,
8,12,6,2,6,12,8,0,
128,128,128,128,128,128,128,128,
0,0,3,7,4,0,0,0,
32,116,84,84,84,124,120,0,
127,127,36,68,68,124,56,0,
56,124,68,68,68,76,72,0,
56,124,68,68,36,127,127,0,
56,124,84,84,84,92,88,0,
0,4,126,127,5,5,1,0,
24,188+128,164+128,164+128,148+256,252,124,0, /* g */
127,127,4,4,4,124,120,0,
0,0,68,125,125,64,0,0,
128,256,132+128,253+256,125+128,0,0,0, /* j */
127,127,16,56,108,68,0,0,
0,0,65,127,127,64,0,0,
124,124,8,24,12,124,120,0,
124,124,4,4,4,124,120,0,
56,124,68,68,68,124,56,0,
252+256,252+256,36,68,68,124,56,0, /* p */
56,124,68,68,36,252+256,252+256,0, /* q */
0,124,124,8,4,4,4,0,
72,92,84,84,84,116,32,0,
0,4,63,127,68,68,0,0,
60,124,64,64,64,124,124,0,
28,60,96,64,96,60,28,0,
60,124,96,48,96,124,60,0,
68,108,56,16,56,108,68,0,
28,188+128,160+128,160+128,144+256,252,124,0, /* y */
68,100,116,92,76,68,0,0,
0,8,8,62,119,65,65,0,
0,0,119,119,0,0,0,0,
0,65,65,119,62,8,8,0,
4,6,2,6,4,6,2,0,
120,124,70,67,70,124,120,0,
68,206,155,177,241,91,10,0,
58,122,64,64,32,122,122,0,
56,124,84,84,85,93,25,0,
34,119,85,85,125,123,66,0,
33,117,84,84,124,121,65,0,
33,117,85,84,124,120,64,0,
32,116,87,87,127,120,64,0,
8,92,212,180,244,84,0,0,
58,127,85,85,85,95,26,0,
57,125,84,84,84,93,25,0,
57,125,85,84,84,92,24,0,
0,1,69,124,124,65,1,0,
2,3,69,125,125,67,2,0,
1,1,69,124,124,64,0,0,
121,125,22,18,22,125,121,0,
112,120,43,43,43,120,112,0,
124,124,84,84,85,69,69,0,
40,116,84,56,124,84,88,0,
126,127,9,127,127,73,73,0,
50,123,73,73,73,123,50,0,
50,122,72,72,72,122,50,0,
50,122,74,72,72,120,48,0,
58,123,65,65,33,123,122,0,
58,122,66,64,32,120,120,0,
26,186,160,160,144,250,122,0,
25,61,102,66,102,61,25,0,
61,125,64,64,64,125,61,0,
24,60,36,126,126,36,36,0,
40,126,127,41,67,114,48,0,
1,11,46,124,124,46,11,1,
127,127,5,21,125,127,82,0,
64,200,136,254,127,11,10,0,
32,116,84,84,125,121,65,0,
0,0,68,125,125,65,0,0,
48,120,72,72,74,122,50,0,
56,120,64,96,58,122,66,0,
74,114,122,10,10,122,112,0,
125,125,25,17,33,125,125,0,
0,18,23,21,23,22,20,0,
2,23,21,21,21,23,2,0,
0,48,120,77,69,96,32,0,
0,56,56,8,8,8,0,0,
0,8,8,8,56,56,0,0,
7,23,184,188,214,242,160,0,
7,87,120,108,246,242,64,0,
0,0,48,125,125,48,0,0,
8,28,54,42,28,54,34,0,
34,54,28,42,54,28,8,0,
170,0,85,0,170,0,85,0,
170,85,170,85,170,85,170,85,
85,255,170,255,85,255,170,255,
0,0,0,255,255,0,0,0,
16,16,16,255,255,0,0,0,
20,20,20,255,255,0,0,0,
16,16,255,255,0,255,255,0,
16,16,240,240,16,240,240,0,
20,20,20,252,252,0,0,0,
20,20,247,247,0,255,255,0,
0,0,255,255,0,255,255,0,
20,20,244,244,4,252,252,0,
20,20,23,23,16,31,31,0,
16,16,31,31,16,31,31,0,
20,20,20,31,31,0,0,0,
16,16,16,240,240,0,0,0,
0,0,0,31,31,16,16,16,
16,16,16,31,31,16,16,16,
16,16,16,240,240,16,16,16,
0,0,0,255,255,16,16,16,
16,16,16,16,16,16,16,16,
16,16,16,255,255,16,16,16,
0,0,0,255,255,20,20,20,
0,0,255,255,0,255,255,16,
0,0,31,31,16,23,23,20,
0,0,252,252,4,244,244,20,
20,20,23,23,16,23,23,20,
20,20,244,244,4,244,244,20,
0,0,255,255,0,247,247,20,
20,20,20,20,20,20,20,20,
20,20,247,247,0,247,247,20,
20,20,20,23,23,20,20,20,
16,16,31,31,16,31,31,16,
20,20,20,244,244,20,20,20,
16,16,240,240,16,240,240,16,
0,0,31,31,16,31,31,16,
0,0,0,31,31,20,20,20,
0,0,0,252,252,20,20,20,
0,0,240,240,16,240,240,16,
16,16,255,255,16,255,255,16,
20,20,20,255,255,20,20,20,
16,16,16,31,31,0,0,0,
0,0,0,240,240,16,16,16,
255,255,255,255,255,255,255,255,
240,240,240,240,240,240,240,240,
255,255,255,255,0,0,0,0,
0,0,0,0,255,255,255,255,
15,15,15,15,15,15,15,15,
56,124,68,56,56,108,68,0,
126,127,21,37,47,58,16,0,
66,126,126,2,2,2,6,0,
2,126,126,2,126,126,2,0,
99,119,93,73,65,99,99,0,
60,126,66,102,62,26,2,0,
64,126,62,32,32,62,30,0,
4,6,66,126,124,6,2,0,
73,93,119,99,119,93,73,0,
28,62,107,73,107,62,28,0,
76,126,115,1,115,126,76,0,
48,122,79,69,101,57,25,0,
24,60,36,60,60,36,60,24,
92,126,58,46,38,63,29,0,
28,62,107,73,73,0,0,0,
126,127,1,1,1,127,126,0,
42,42,42,42,42,42,42,0,
0,68,68,95,95,68,68,0,
0,64,81,91,78,68,64,0,
0,64,68,78,91,81,64,0,
0,0,0,254,255,3,2,0,
0,32,96,127,63,0,0,0,
0,16,16,84,84,16,16,0,
36,54,18,54,36,54,18,0,
6,15,9,9,9,15,6,0,
0,0,0,24,24,0,0,0,
0,0,0,16,16,0,0,0,
16,16,48,127,127,1,1,1,
1,15,14,1,15,14,0,0,
10,11,13,15,10,0,0,0,
0,60,60,60,60,60,0,0,
0,0,0,0,0,0,0,0};