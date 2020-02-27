#----------------- LIBRARY -------------------
BINDIR = /home/skameno/Programs/PolaR/
#----------------- LINK OPTIONS -------------------
CCOMPL=gcc $(CFLAGS)
#------- Followings are PASS or DIRECTORY -------
PROGS=	PolariSplit PolariBunch PolariTimeShift SpecInteg
GRLIBS= -L/usr/X11R6/lib -lX11
MATH=	-lm
#----------------- MAPPING ------------------------
OBJ_SpecInteg  =	SpecInteg.o
OBJ_PolariSplit=	PolariSplit.o
OBJ_PolariBunch=	PolariBunch.o
OBJ_TimeShift  =	PolariTimeShift.o
#----------------- Compile and link ------------------------
PolariSplit : $(OBJ_PolariSplit)
	$(CCOMPL) -o $@ $(OBJ_PolariSplit)

PolariBunch : $(OBJ_PolariBunch)
	$(CCOMPL) -o $@ $(OBJ_PolariBunch)

PolariTimeShift : $(OBJ_TimeShift)
	$(CCOMPL) -o $@ $(OBJ_TimeShift)

SpecInteg : $(OBJ_SpecInteg)
	$(CCOMPL) -o $@ $(OBJ_SpecInteg)

clean :
	\rm $(PROGS) *.o a.out core *.trace

all :	$(PROGS)

install:
	@mv $(PROGS) $(BINDIR)

#----------------- Objects ------------------------
.c.o:
	$(CCOMPL) -c $*.c
PolariSplit.o:      PolariSplit.c       shm_k5data.inc
PolariBunch.o:	PolariBunch.c shm_k5data.inc
SpecInteg.o:	SpecInteg.c shm_k5data.inc
TimeShift.o:	PolariTimeShift.c shm_k5data.inc
#----------------- End of File --------------------
