LIB = libRES_Event.so
DICT = RES_EventDict

LD = g++

CXXFLAGS += -fPIC -I. -I$(ROOTSYS)/include

all: $(LIB)

clean:
	@echo "Cleaning up RES_Event ..."
	@rm -f *.o $(DICT).* $(LIB)

$(LIB): $(DICT).o RES_Event.o
	$(LD) -o $(LIB) $(DICT).o RES_Event.o -shared
	@echo "DONE."

$(DICT).cc: RES_Event.hh RES_Event.cc
	rootcint -f $(DICT).cc -c -p RES_Event.hh LinkDef.h