# ***************************************************************************
  CC = gcc 
  FF = gcc
  STD = ANSI
  ARCH = Sparc
  # CFLAGS = -O 
  CFLAGS = -O3 
  FLAGS  = -O3
  DataType = BIN 
# ***************************************************************************
#
# There are two paths that you need to change in this Makefile.
# The first is the Genitor root directory and the second is the
# Lcl directory about 23 lines later.

# Genitor root directory 
#GaRoot = /home/oshir/CHC_CSUga/
GaRoot  = /home/naamah/Documents/CHC/CHC_CSUga/

# Fortran Linkage
LIBS = -lgfortran 

# Library Directories
  GaLibDir     = $(GaRoot)Lib/$(ARCH)/

# Actual Libraries to be Used
  GaLib    = Gacsu_$(DataType)

# Include Directories
  GaInc = $(GaRoot)Ga/Include/
  SptInc = $(GaRoot)CHC/Include/

# Src Support Directories
  SptSrc = $(GaRoot)CHC/Src/

# evaluation code directories
  EvSrc = $(GaRoot)CHC/Evals/
  EvInc = $(GaRoot)CHC/Evals/

# Place where your executable will be built 
  #Lcl = $(GaRoot)CHC/
  Lcl = $(GaRoot)CHC/

# *****************************************
# application specific, data type dependent 
# *****************************************
MainProgram = main
      EvMod = m_evals
        SHG = shg
    OpName1 = op_hux
    OpName2 = op_cataclysm
   ProgName = XCHC
   BinSpt   = binarys
   StatSpt  = stats

FOBJS = numrecipes.o shaper.o sleep.o ch.o cclock.o

# ********************************
# let's make an executable program
# ********************************
all: $(ProgName)

$(ProgName) : $(MainProgram).o $(EvMod).o $(FOBJS) $(SHG).o $(OpName1).o $(OpName2).o \
	      $(BinSpt).o $(StatSpt).o
	$(CC) $(CFLAGS) -o $(ProgName)  $(MainProgram).o $(EvMod).o $(FOBJS) $(SHG).o \
	      $(OpName1).o $(OpName2).o $(BinSpt).o $(StatSpt).o \
	      -L$(GaLibDir) -l$(GaLib) $(LIBS) -lm 
	rm -rf $(Lcl)*.o
  
# *******************************************************
# application independent, data type independent function
# information to include in main() driver file
# *******************************************************
GaHeaders = $(GaInc)ga_random.h  $(GaInc)gene.h      $(GaInc)ga_status.h \
	    $(GaInc)ga_params.h  $(GaInc)ga_pool.h   $(GaInc)ga_selection.h \
	    $(GaInc)ga_signals.h $(GaInc)ga_global.h $(GaInc)ga_xdr.h 

# *****************************************
# application specific main driver function
# *****************************************
$(MainProgram).o : $(MainProgram).c $(GaHeaders) $(EvInc)$(EvMod).h \
		   $(SptInc)binarys_params.h $(EvInc)dejong_params.h \
		   $(SptInc)$(OpName1).h $(SptInc)$(OpName2).h  \
		   $(SptInc)op_hux_params.h
	$(CC) -c $(CFLAGS) $(MainProgram).c -I$(GaInc) -I$(SptInc) -I$(EvInc) \
	      -D$(DataType) -D$(STD) -o $(Lcl)$(MainProgram).o

# ***********************************************************
# application specific evaluation function 
# ***********************************************************
$(EvMod).o : $(EvSrc)$(EvMod).c $(GaInc)gene.h $(EvInc)$(EvMod).h \
	     $(EvInc)dejong_params.h
	$(CC) -c $(CFLAGS) $(EvSrc)$(EvMod).c -I$(GaInc) -I$(EvInc) \
	      -I$(SptInc) -D$(DataType) -D$(STD) -o $(Lcl)$(EvMod).o

$(SHG).o : $(EvSrc)$(SHG).F 
	$(FF) -c $(FFLAGS) $(EvSrc)$(SHG).F -I$(EvInc) -o $(Lcl)$(SHG).o

numrecipes.o: $(EvSrc)numrecipes.f
	$(FF) $(FFLAGS) -I./ -c $<
shaper.o: $(EvSrc)shaper.f
	$(FF) $(FFLAGS) -I./ -c $<
ch.o: $(EvSrc)ch.f
	$(FF) $(FFLAGS) -I./ -c $<
cclock.o: $(EvSrc)cclock.c
	$(CC) $(CFLAGS) -I./ -c $<
sleep.o: $(EvSrc)sleep.c
	$(CC) $(CFLAGS) -I./ -c $<

# ***********************************************************
# HUX Operator
# ***********************************************************
$(OpName1).o : $(SptSrc)$(OpName1).c $(GaInc)gene.h $(SptInc)$(OpName1).h 
	$(CC) -c $(CFLAGS) $(SptSrc)$(OpName1).c -I$(GaInc) \
	      -I$(SptInc) -D$(DataType) -D$(STD) -o $(Lcl)$(OpName1).o


# ***********************************************************
# Cataclysmic Mutation Operator
# ***********************************************************
$(OpName2).o : $(SptSrc)$(OpName2).c $(GaInc)gene.h $(SptInc)$(OpName2).h 
	$(CC) -c $(CFLAGS) $(SptSrc)$(OpName2).c -I$(GaInc) \
	      -I$(SptInc) -D$(DataType) -D$(STD) -o $(Lcl)$(OpName2).o


# ***********************************************************
# Binary Application Support
# ***********************************************************
$(BinSpt).o : $(SptSrc)$(BinSpt).c $(GaInc)gene.h $(SptInc)$(BinSpt).h  \
	      $(SptInc)binarys_params.h
	$(CC) -c $(CFLAGS) $(SptSrc)$(BinSpt).c -I$(GaInc) \
	      -I$(SptInc) -D$(DataType) -D$(STD) -o $(Lcl)$(BinSpt).o


# ***********************************************************
# Statistical Support
# ***********************************************************
$(StatSpt).o : $(SptSrc)$(StatSpt).c $(GaInc)gene.h $(SptInc)$(StatSpt).h  \
	       $(SptInc)stats_params.h
	$(CC) -c $(CFLAGS) $(SptSrc)$(StatSpt).c -I$(GaInc) \
	      -I$(SptInc) -D$(DataType) -D$(STD) -o $(Lcl)$(StatSpt).o

