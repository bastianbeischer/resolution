name := resolution
G4TARGET := $(name)
G4EXLIB := true

.PHONY: all
all: lib bin TAGS

CPPFLAGS += -I$(ROOTSYS)/include
EXTRALIBS += $(shell root-config --libs)

include $(G4INSTALL)/config/binmake.gmk

TAGS: $(G4TARGET).cc src/*.cc include/*.hh
	@rm -f TAGS;
	etags $(G4TARGET).cc src/*.cc include/*.hh
