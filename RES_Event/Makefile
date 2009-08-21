LIB = libRES_Event.so
DICT = RES_EventDict

LD = g++

CXXFLAGS	+= -I. -I$(ROOTSYS)/include

all: $(LIB)

clean:
		@rm -f *.o $(DICT).* $(LIB)



$(LIB): $(DICT).o RES_Event.o
		@echo "Creating shared lib libRES_Event.so"
		$(LD) -o $(LIB) $(DICT).o RES_Event.o -shared
		@echo "DONE."
		@echo

$(DICT).cc: RES_Event.hh RES_Event.cc
		rootcint -f $(DICT).cc -c -p RES_Event.hh LinkDef.h