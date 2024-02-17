static char Pgeo[]="\n\
Environment variables SHOWGEOMETRY, PLOTGEOMETRY, BLENDGEOMETRY, JKVGEOMETRY,\n\
and the arguments of -g (warning: the syntax for blend is modified) accept:\n\
Standard X-geometry format:\n\
  XSIZExYSIZE = define size of the window (without decorations) in pixels\n\
  XSIZExYSIZE[+-]DX[+-]DY = as above with position\n\
    +=top/left, -=bottom/right (bug: does not handle correctly decoration sizes)\n\
Extensions (keywords can be shortened to 1 letter):\n\
  SIZE = SIZExSIZE\n\
  f = fullscreen = BIGxBIG = full screen mode\n\
  v = vga = 640x480    g = 800x600    b = beamer = 1024x768\n\
  n = 1680x1050\n\
  h = high-definition = 1920x1080     r = 1280×720 (HD ready)\n\
  p = portrait = 768x1024             P = Portrait = 960x1050\n\
  s = square = 600x600+0+0            S = Square = 1024x1024\n\
   = use the size of the first picture (jkv only: jkv -g pic.png)\n\
This applies to Ubuntu Mate:\n\
  A window can be enlarged, not shrank (to be reconsidered..).\n\
  Alt-F10 fullscreen is without the area occupied by the panel (incl. hidden)\n\
    and titlebar (unless \"Undecorate maximized widows\" in MATE Tweak).\n\
  True fullscreen (-gf) blocks the panel (but not further raised windows).\n\
";
