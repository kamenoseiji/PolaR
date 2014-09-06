#----------------- LIBRARY -------------------
BINDIR = /usr/custom/bin
#----------------- LINK OPTIONS -------------------
CCOMPL=gcc $(CFLAGS)
#------- Followings are PASS or DIRECTORY -------
PROGS=	PolariSplit PolariBunch
GRLIBS= -L/usr/X11R6/lib -lX11
MATH=	-lm
#----------------- MAPPING ------------------------
OBJ_PolariBunch=	PolariBunch.o
#----------------- Compile and link ------------------------
PolariBunch : $(OBJ_PolariBunch)
	$(CCOMPL) -o $@ $(OBJ_PolariBunch)

clean :
	\rm $(PROGS) *.o a.out core *.trace

all :	$(PROGS)

install:
	@mv $(PROGS) $(BINDIR)

#----------------- Objects ------------------------
.c.o:
	$(CCOMPL) -c $*.c
PolariBunch.o:	PolariBunch.c shm_k5data.inc
#----------------- End of File --------------------
