LIB = libAMS02_Magnet.so

LD = g++

CXXFLAGS += -fPIC -I. -I$(ROOTSYS)/include -D_PGTRACK_

all: $(LIB)

clean:
	@echo "Cleaning up AMS02_Magnet ..."
	@rm -f *.o $(LIB)

$(LIB): MagField.o
	$(LD) -o $(LIB) MagField.o -shared `root-config --libs`
