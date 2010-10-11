STEM = RES_Event

TARGET = lib$(STEM).so
DICT = $(STEM)Dict

LD = g++

CXXFLAGS += -fPIC -I. -I$(shell root-config --incdir)
LIBS += `root-config --libs`

all: $(TARGET)

clean:
	@echo "Cleaning up $(STEM) ..."
	@rm -f *.o $(DICT).* $(TARGET)

$(TARGET): $(DICT).o $(STEM).o
	$(LD) -o $(TARGET) $(DICT).o $(STEM).o -shared $(LIBS)

$(DICT).cc: $(STEM).hh $(STEM).cc
	rootcint -f $(DICT).cc -c -p $(STEM).hh LinkDef.h
