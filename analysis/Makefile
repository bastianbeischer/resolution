CC = g++

LDFLAGS += $(shell root-config --libs) -L../RES_Event -lRES_Event
CXXFLAGS += -I$(ROOTSYS)/include -I../RES_Event -I$(MY_INCLUDES)

.PHONY: all clean

all: single_run res_vs_mom res_vs_angle chi2_dist

single_run: single_run.o
res_vs_mom: res_vs_mom.o
res_vs_angle: res_vs_angle.o
chi2_dist: chi2_dist.o

clean:
	@echo "Cleaning up analysis ..."
	@rm -f *.o single_run res_vs_mom res_vs_angle chi2_dist