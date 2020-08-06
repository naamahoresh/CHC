# ***************************************************************************
# 
# This makefile allows you to create libraries of fundamental Genitor 
# functions (files containing these functions begin with `ga_'.  
# You should make libraries for each gene type to be used, for each 
# architecture that it is used on, and if prototyping is to be used or not
# (i.e., ANSI or NON_ANSI respectively).
#
# The libraries can be built by the command:
#
# make clean STD=ANSI CC=gcc CFLAGS=-O DATA=INT ARCH=Sparc
#
# where the valid options are:
#
#        STD = ANSI or NON_ANSI
#        CC = cc or gcc (other compilers have not been tested but may be used)
#        CFLAGS = -O or -g and (-I/usr/ucbinclude for Solaris)
#        DATA = INT, BIN, FLOAT
#        ARCH = Solaris, Sparc, Aix, Dec, HP
#
#
# Although your applications must link with the appropriate library for the
# data type to be used, the function calls themselves should be identical
# regardless of the data type.  See the GENITOR CODE STRUCTURE document
# for details.
#
# ***************************************************************************
STD=ANSI
CC=gcc
CFLAGS=-O3 -ansi
DATA=BIN
ARCH=Sparc

#**************************
# where are things located?
# *************************

# GaRoot is the only thing that you should want to alter.
#GaRoot  = /home/oshir/CHC_CSUga/
GaRoot  = /home/naamah/Documents/CHC/CHC_CSUga/

GaSrc   = $(GaRoot)Ga/Src/
GaInc   = $(GaRoot)Ga/Include/
GaLib   = $(GaRoot)Lib/$(ARCH)/


# ********************************************************
# application independent, data type independent functions
# ********************************************************
GaObjects = $(GaSrc)ga_tournament.o $(GaSrc)ga_xdr.o $(GaSrc)ga_params.o \
	    $(GaSrc)ga_signals.o $(GaSrc)ga_status.o $(GaSrc)ga_pool.o \
	    $(GaSrc)ga_selection.o $(GaSrc)ga_copy.o $(GaSrc)ga_random.o

# ************************
# library of ga_ functions
# ************************

LibName = $(GaLib)libGacsu_$(DATA)

library : $(LibName).a

$(LibName).a: $(GaObjects)
	@if test $(ARCH) = "Solaris" ;                      \
	then                                                \
	  /usr/ccs/bin/ar rv $(LibName).a $(GaObjects);     \
	  /bin/true $(LibName).a;                           \
	else                                                \
	  ar rv $(LibName).a $(GaObjects);                  \
	  ranlib $(LibName).a;                              \
	fi
	rm -rf $(GaSrc)*.o

clean:
	rm -rf $(GaSrc)*.o
	rm -rf $(LibName).a 
	make 



# ********************************************************
# application independent, data type independent functions
# ********************************************************
$(GaSrc)ga_tournament.o : $(GaSrc)ga_tournament.c $(GaInc)ga_tournament.h \
			  $(GaInc)gene.h  
	$(CC) -c $(CFLAGS) $(GaSrc)ga_tournament.c -D$(DATA) -D$(STD) \
	   -I$(GaInc) -o $(GaSrc)ga_tournament.o

$(GaSrc)ga_xdr.o : $(GaSrc)ga_xdr.c $(GaInc)ga_xdr.h $(GaInc)gene.h 
	$(CC) -c $(CFLAGS) $(GaSrc)ga_xdr.c -D$(DATA) -D$(STD) \
	   -I$(GaInc) -o $(GaSrc)ga_xdr.o

$(GaSrc)ga_params.o : $(GaSrc)ga_params.c $(GaInc)ga_params.h $(GaInc)gene.h \
		      $(GaInc)ga_global_extern.h $(GaInc)ga_status.h 
	$(CC) -c $(CFLAGS) $(GaSrc)ga_params.c -D$(DATA) -D$(STD) \
	   -I$(GaInc) -o $(GaSrc)ga_params.o

$(GaSrc)ga_signals.o : $(GaSrc)ga_signals.c $(GaInc)ga_status.h 
	$(CC) -c $(CFLAGS) $(GaSrc)ga_signals.c -D$(DATA) -D$(STD) \
	   -I$(GaInc) -o $(GaSrc)ga_signals.o

$(GaSrc)ga_status.o : $(GaSrc)ga_status.c $(GaInc)ga_status.h $(GaInc)gene.h \
		      $(GaInc)ga_global_extern.h $(GaInc)ga_pool.h 
	$(CC) -c $(CFLAGS) $(GaSrc)ga_status.c -D$(DATA) -D$(STD) \
	   -I$(GaInc) -o $(GaSrc)ga_status.o

$(GaSrc)ga_pool.o : $(GaSrc)ga_pool.c $(GaInc)ga_status.h $(GaInc)gene.h \
		    $(GaInc)ga_random.h $(GaInc)ga_pool.h $(GaInc)ga_copy.h 
	$(CC) -c $(CFLAGS) $(GaSrc)ga_pool.c -D$(DATA) -D$(STD) \
	   -I$(GaInc) -o $(GaSrc)ga_pool.o

$(GaSrc)ga_selection.o : $(GaSrc)ga_selection.c $(GaInc)gene.h \
			 $(GaInc)ga_random.h $(GaInc)ga_copy.h 
	$(CC) -c $(CFLAGS) $(GaSrc)ga_selection.c -D$(DATA) -D$(STD) \
	   -I$(GaInc) -o $(GaSrc)ga_selection.o

$(GaSrc)ga_copy.o : $(GaSrc)ga_copy.c $(GaInc)ga_copy.h $(GaInc)gene.h 
	$(CC) -c $(CFLAGS) $(GaSrc)ga_copy.c -D$(DATA) -D$(STD) \
	   -I$(GaInc) -o $(GaSrc)ga_copy.o

$(GaSrc)ga_random.o : $(GaSrc)ga_random.c $(GaInc)ga_random.h $(GaInc)gene.h 
	$(CC) -c $(CFLAGS) $(GaSrc)ga_random.c -D$(DATA) -D$(STD) \
	   -I$(GaInc) -o $(GaSrc)ga_random.o
