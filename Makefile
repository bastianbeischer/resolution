name := resolution
G4TARGET := $(name)
G4EXLIB := true

CPPFLAGS += -I$(ROOTSYS)/include
EXTRALIBS += $(shell root-config --libs)

RES_EventDir := ./RES_Event
RES_EventLib := $(RES_EventDir)/libRES_Event.so
EXTRALIBS += -L./$(RES_EventDir) -lRES_Event
CPPFLAGS += -I$(RES_EventDir)
EXTRA_LINK_DEPENDENCIES += $(RES_EventLib)

BlobelDir := ./blobel
BlobelLib := $(BlobelDir)/libBlobel.so
EXTRALIBS += -L./$(BlobelDir) -lBlobel
CPPFLAGS += -I$(BlobelDir)
EXTRA_LINK_DEPENDENCIES += $(BlobelLib)

.PHONY: all
all: lib bin TAGS

include $(G4INSTALL)/config/binmake.gmk

TAGS: $(G4TARGET).cc src/*.cc include/*.hh
	@rm -f TAGS;
	etags $(G4TARGET).cc src/*.cc include/*.hh

$(RES_EventLib): $(RES_EventDir)/*.hh $(RES_EventDir)/*.cc
	$(MAKE) -C $(RES_EventDir)

clean_ams_event:
	$(MAKE) -C $(RES_EventDir) clean

$(BlobelLib): $(BlobelDir)/*.h $(BlobelDir)/*.f 
	$(MAKE) -C $(BlobelDir)

test: test.cc
	g++ -o test.o -I$(ROOTSYS)/include -IRES_Event -c test.cc
	g++ -o test $(shell root-config --libs) -LRES_Event -l RES_Event test.o

clean_test:
	rm -f test.o test