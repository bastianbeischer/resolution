# Some directories
SOURCEDIR    := src
INCLUDEDIR   := include
TEMPDIR      := .tmp

# Name of the output executable
EXE          := analysis

# Specify the input files
MAIN         := main
SOURCES      := SingleFile.cc ResVsMom.cc
HEADERS      := SingleFile.hh ResVsMom.hh
OBJECTS      := $(TEMPDIR)/$(MAIN).o $(addprefix $(TEMPDIR)/,$(SOURCES:.cc=.o))

# Compilers and Flags
CC           := gcc
CXX          := g++
INCLUDEFLAGS := -I$(INCLUDEDIR) -I$(shell root-config --incdir) -I$(MY_INCLUDES) -I../RES_Event
LIBS         := $(shell root-config --libs) -L../RES_Event -lRES_Event
CCFLAGS      := -Wall -fPIC $(INCLUDEFLAGS) -MMD
CXXFLAGS     := -Wall -fPIC $(INCLUDEFLAGS) -MMD

# Documenation
DOCDIR       := $(PWD)/doc
DOXYGEN_FILE := $(EXE).doxy

# What do we want to build?
all: exe TAGS

# PHONY targets
.PHONY: all clean dist_clean doc exe check-syntax

clean:
	@echo Cleaning up...
	@rm -fr $(TEMPDIR)

dist_clean:
	@echo Cleaning up everything...
	@rm -fr $(TEMPDIR) $(DOCDIR) $(EXE)  TAGS

doc:
	@echo Creating documentation
	@doxygen $(DOXYGEN_FILE)	

exe: $(EXE)

check-syntax:
	$(CXX) $(CXXFLAGS) -fsyntax-only $(CHK_SOURCES)

# Regular build targets
$(EXE): $(TEMPDIR) $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJECTS) $(LIBS)

$(TEMPDIR):
	@if [ ! -d $(TEMPDIR) ]; then mkdir $(TEMPDIR); fi

$(TEMPDIR)/%.o: $(SOURCEDIR)/%.cc
	$(CXX) $(CXXFLAGS) -o $@ -c $<

$(TEMPDIR)/%.o: $(SOURCEDIR)/%.c
	$(CC) $(CCFLAGS) -o $@ -c $<

$(TEMPDIR)/%.o: %.cc
	$(CXX) $(CXXFLAGS) -o $@ -c $<

TAGS: $(MAIN).cc $(addprefix $(SOURCEDIR)/,$(SOURCES)) $(addprefix $(INCLUDEDIR)/, $(HEADERS))
	@rm -f TAGS;
	etags $^

# Include automatically regenerated dependency files, if they exist
-include $(OBJECTS:.o=.d)
