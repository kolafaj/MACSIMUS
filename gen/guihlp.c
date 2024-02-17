static char Pgui[]="\n\
Environment variable GUI (common setup of the graphical interface):\n\
  GUI = \"KEYS dDELAY #HEXCOL\"  (in this order, spaces optional)\n\
where\n\
  KEYS = String containg the following chars:    +---------------------------+\n\
    s = `show' will be started with menu on      + F10 or shaded menu button +\n\
    p = `plot' will be started with menu on      +    toggles menu on/off    +\n\
    b = `blend' will be started with menu on     +---------------------------+\n\
    K = `plot' will enable button [KILL ALL] (hotkey K) calling `killall plot'\n\
        Does not apply for button [kill all] (hotkey K) available in plots\n\
        spawn by `showcp' or `rdfg' which kills the spawn plots only.\n\
    v = verbose printout to stderr\n\
  dDELAY = Delay in s before redraw for (a series of) ConfigureNotify\n\
           or Expose events, range=0.001-2, default=0.2.\n\
           Increase if the connection is slow.\n\
  #HEXCOL = Background (\"black\") color (in hexadecimal web style).\n\
            Only dark colors are recommended.\n\
            Does not apply to `show' in 3D modes (cf. `show -bg')\n\
Example:\n\
  export GUI=\"bsK d0.5 #003010\"\n\
";
