CC = g++

LDFLAGS += $(shell root-config --libs) -L../RES_Event -lRES_Event
CXXFLAGS += -I$(ROOTSYS)/include -I../RES_Event -I$(MY_INCLUDES)

TARGETS = single_run res_vs_mom res_vs_angle chi2_dist compare_msc

.PHONY: all clean

all: $(TARGETS)

%: %.o

clean:
	@echo "Cleaning up analysis ..."
	@rm -f *.o $(TARGETS)