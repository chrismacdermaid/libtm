.PHONY: clean clean-all install install-latest archive update

VMDPLUGINS = $(addprefix $(HOME)/,.vmdplugins)

clean-all:
	cd src; make clean-all 

clean-%:
	cd src; make clean-$(@:clean-%=%);

makeshlib-%:
	cd src; make makeshlib
	cd src; make -f Makefile.shlib $(@:makeshlib-%=%);

archive:
	git archive --prefix=tclmath/ HEAD -o tclmath-latest.zip 
