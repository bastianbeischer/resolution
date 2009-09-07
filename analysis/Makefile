CC = g++
LD = g++

LDFLAGS += $(shell root-config --libs) -L../RES_Event -lRES_Event
CXXFLAGS += -I$(ROOTSYS)/include -I../RES_Event -I$(MY_INCLUDES)

.PHONY: all clean

all: analysis analysis2

analysis: analysis.o
analysis2: analysis2.o

clean:
	rm -f *.o analysis analysis2