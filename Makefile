# $Id: Makefile,v 1.3 2004/12/02 17:04:01 bittner Exp $
# $Revision: 1.3 $
# $Date: 2004/12/02 17:04:01 $
# $Author: bittner $
#Makefile,v
#Revision 1.2  2004/09/29 00:48:59  bittner
#
#Added altvec compiler options for OSX machines.
#
#Revision 1.1  2004/09/24 23:08:21  bittner
#Initial revision
#
#  MAKE FILE FOR EXCITON CODE

.SUFFIXES: .f .F


EXCITON_BIN_DIR = ~/bin

OBJ = eigvals.o main-exci.o input.o integrate.o io.o  \
      prepare.o propers.o relaxation.o tools.o

INCLUDES = constant.fi coords.fi elstates.fi parameters.fi sizes.fi\
            structure.fi couplings.fi genparam.fi rates.fi spectra.fi  vibmodes.fi

CODE = eigvals.f main-exci.f input.f integrate.f io.f  \
      prepare.f propers.f relaxation.f tools.f


# xlf compiler
#FCFLAGS = -c -qfixed=132
#F90  = /opt/ibmcmp/xlf/8.1/bin/xlf90 $(FCFLAGS)
#FC   = /opt/ibmcmp/xlf/8.1/bin/xlff90
#F90FLAGS = -c -w -O -W 132 -YEXT_NAMES=LCS -YEXT_SFX=_  -g -TENV:simd_imask=off 
#LKFLAGS = -X -framework -X vecLib -lU77 
#LKFLAGS = -lU77 

#F90 = /Applications/Absoft13.0/bin/f90 $(F90FLAGS)
#FC = /Applications/Absoft13.0/bin/f90



#------ GNU FORTRAN COMPILER------


FCFLAGS = -c -w  -O  -ffixed-line-length-none -fbounds-check
F90 = /usr/bin/gfortran $(FCFLAGS)
FC = /usr/bin/gfortran 



DISTRIB_LEVEL = 2.0.1
DISTRIB_DIR = exciton-$(DISTRIB_LEVEL)
DISTRIB_SRC = $(CODE) $(INCLUDES) Makefile 
DOC = README INSTALL GPL.txt
DISTRIB_EXAM = $(DISTRIB_DIR)/test
DISTRIB_SRCDIR = $(DISTRIB_DIR)/src

all:
	make exciton

install:
	make all
	cp exciton.x $(EXCITON_BIN_DIR)/bin



distrib: $(DISTRIB_SRC) $(DOC) Makefile test
	mkdir $(DISTRIB_DIR)
	mkdir $(DISTRIB_SRCDIR)
	mkdir $(DISTRIB_EXAM)
	cp $(DISTRIB_SRC) $(DISTRIB_SRCDIR)
	cp -Rf test  $(DISTRIB_DIR)
	cp $(DOC) $(DISTRIB_DIR)
	tar cvf exciton-$(DISTRIB_LEVEL).tar exciton-$(DISTRIB_LEVEL)
	mv $(DISTRIB_DIR).tar distribs
	rm -rf $(DISTRIB_DIR)

electrode.o: electrode.f
	$(F90) electrode.f

eigvals.o: sizes.fi constant.fi eigvals.f
	$(F90) eigvals.f

main-exci.o: $(INCLUDES) main-exci.f
	$(F90) main-exci.f

input.o: sizes.fi structure.fi constant.fi input.f
	$(F90) input.f	

integrate.o: elstates.fi parameters.fi sizes.fi spectra.fi integrate.f
	$(F90) integrate.f

io.o: $(INCLUDES) io.f
	$(F90) io.f

prepare.o: $(INCLUDES) prepare.f
	$(F90) prepare.f

propers.o: $(INCLUDES) propers.f
	$(F90) propers.f

relaxation.o: $(INCLUDES) relaxation.f
	$(F90) relaxation.f

tools.o: $(INCLUDES) tools.f
	$(F90) tools.f

exciton: $(OBJ)
	$(FC) -o exciton.x $(OBJ)  $(LKFLAGS)


clean:
	rm -f *.o
# RULES
.f.o :
	$(F90) $<	
.F.o :
	$(F90) $<
