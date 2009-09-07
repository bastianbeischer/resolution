name := resolution
G4TARGET := $(name)
G4EXLIB := true

#CPPVERBOSE := 1
#CPPFLAGS += -Wno-unused-result

CPPFLAGS += -I$(ROOTSYS)/include
EXTRALIBS += $(shell root-config --libs) -lMinuit -lMinuit2

CPPFLAGS += -g

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

TAGS: $(G4TARGET).cc src/*.cc include/*.hh $(RES_EventDir)/*.cc $(RES_EventDir)/*.hh $(BlobelDir)/*.h $(BlobelDir)/*.f
	@rm -f TAGS;
	etags $^

clean:: clean_res_event clean_blobel clean_analysis clean_test

$(RES_EventLib): $(RES_EventDir)/*.hh $(RES_EventDir)/*.cc
	$(MAKE) -C $(RES_EventDir)

clean_res_event:
	$(MAKE) -C $(RES_EventDir) clean

$(BlobelLib): $(BlobelDir)/*.h $(BlobelDir)/*.f 
	$(MAKE) -C $(BlobelDir)

clean_blobel:
	$(MAKE) -C $(BlobelDir) clean

test:
	$(MAKE) -C test

clean_test:
	$(MAKE) -C test clean

analysis:
	$(MAKE) -C analysis

clean_analysis:
	$(MAKE) -C analysis clean