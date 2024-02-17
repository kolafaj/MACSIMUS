#!/bin/bash
rm *.o *.lo *.xo
makemake linux gcc ; mv makefile makefile.double    # mapped to ctrl-t
makemake linux gcc long ; mv makefile makefile.long # mapped to t
makemake linux gcc high ; mv makefile makefile.high # mapped to T
