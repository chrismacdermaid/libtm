# TCLMATH multiple-machine Makefile

# Definitions

ROOT = tclmath 
VERSION = 1.0
EXE =	$(ROOT)_$@
SRC =	$(wildcard *.cpp)
INC =	$(wildcard *.h)
OBJ = 	$(SRC:.cpp=.o)

# List of all targets

help:
	@echo ''
	@echo 'make clean-all           delete all object files'
	@echo 'make clean-machine       delete object files for one machine'
	@echo ''
	@echo 'make makelib             create Makefile.lib'
	@echo 'make makeshlib           create Makefile.shlib'
	@echo 'make -f Makefile.lib     compile static lib: libccb.a'
	@echo 'make -f Makefile.shlib   compile shared object lib: libtm.so required for TCL'     
	@echo '' 
	@echo 'make machine             build tclmath where machine is one of:'
	@echo ''
	@files="`ls MAKE/Makefile.*`"; \
	  for file in $$files; do head -1 $$file; done
	@echo ''

# Build the code

.DEFAULT:
	@test -f MAKE/Makefile.$@
	@if [ ! -d Obj_$@ ]; then mkdir Obj_$@; fi
	@cp -p *.cpp *.h Obj_$@
	@cp MAKE/Makefile.$@ Obj_$@/Makefile
	@cd Obj_$@; \
	$(MAKE) $(MFLAGS) "OBJ = $(OBJ)" "INC = $(INC)" \
	  "EXE = ../$(EXE)" ../$(EXE)
	@if [ -d Obj_$@ ]; then cd Obj_$@; rm -f $(SRC) $(INC) Makefile*; fi

# Remove machine-specific object files
clean:
	@echo 'make clean-all           delete all object files and documentation'
	@echo 'make clean-machine       delete object files for one machine'
	@echo 'make clean-docs          delete documentation files'

clean-all:
	rm -rf Obj_*

clean-%:
	rm -rf Obj_$(@:clean-%=%)

# Update Makefile.lib and Makefile.list
makelib:                                                                                                                                         
	@$(SHELL) Make.sh Makefile.lib


makeshlib:
	@$(SHELL) Make.sh Makefile.shlib

#Emacs on-the-fly syntax checker
check-syntax:
	$(CC) -o /dev/null $(CCB_INC) $(CCFLAGS) -S ${CHK_SOURCES}
