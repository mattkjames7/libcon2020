

ifndef BUILDDIR 
	export BUILDDIR=$(shell pwd)/build
endif

ifeq ($(OS),Windows_NT)
#windows stuff here
	MD=mkdir
else
#linux and mac here
	MD=mkdir -p
endif

.PHONY: all lib obj clean header test

all: obj lib header

windows: winobj winlib header

obj:
	$(MD) $(BUILDDIR)
	cd src; make obj

lib:
	$(MD) lib/libcon2020
	cd src; make lib

winobj:
	$(MD) $(BUILDDIR)
	cd src; make winobj

winlib: 
	$(MD) lib/libcon2020
	cd src; make winlib

header:
ifneq (,$(shell which python3))
	python3 generateheader.py
else
	@echo "python3 command doesn't appear to exist - skipping header regeneration..."
endif

test:
	cd test; make all

clean:
	cd test; make clean
	-rm -v lib/libcon2020/libcon2020.so
	-rm -v lib/libcon2020/libcon2020.dll
	-rmdir -v lib/libcon2020
	-rm -v build/*.o
	-rmdir -v build