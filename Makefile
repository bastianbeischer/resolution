name := resolution
G4TARGET := $(name)
G4EXLIB := true

#CPPVERBOSE := 1

CPPFLAGS += -I$(ROOTSYS)/include
EXTRALIBS += $(shell root-config --libs) -lMinuit

#CPPFLAGS += -g
CPPFLAGS += -Wno-unused-result

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

MillepedeDir := ./millepede
MillepedeLib := $(MillepedeDir)/libMillepede.so
EXTRALIBS += -L./$(MillepedeDir) -lMillepede
CPPFLAGS += -I$(MillepedeDir)
EXTRA_LINK_DEPENDENCIES += $(MillepedeLib)


.PHONY: all
all: TAGS lib bin

include $(G4INSTALL)/config/binmake.gmk

TAGS: $(G4TARGET).cc src/*.cc include/*.hh $(RES_EventDir)/*.cc $(RES_EventDir)/*.hh $(BlobelDir)/*.h $(BlobelDir)/*.f
	@rm -f TAGS;
	@etags $^

clean::
	@$(MAKE) -C $(RES_EventDir) clean
	@$(MAKE) -C $(BlobelDir) clean
	@$(MAKE) -C $(MillepedeDir) clean
	@$(MAKE) -C analysis clean
	@$(MAKE) -C test clean
	@rm -f TAGS;

$(RES_EventLib): $(RES_EventDir)/*.hh $(RES_EventDir)/*.cc
	@echo "Creating shared library $(RES_EventDir)/libRES_Event.so ..."
	@$(MAKE) -s -C $(RES_EventDir)

$(BlobelLib): $(BlobelDir)/*.h $(BlobelDir)/*.f 
	@echo "Creating shared library $(BlobelDir)/libBlobel.so ..."
	@$(MAKE) -s -C $(BlobelDir)

$(MillepedeLib): $(MillepedeDir)/*.h $(MillepedeDir)/*.f 
	@echo "Creating shared library $(MillepedeDir)/libMillepede.so ..."
	@$(MAKE) -s -C $(MillepedeDir)
