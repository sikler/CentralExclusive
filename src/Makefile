ONLXPLUS := $(shell uname -n | grep lx)

ifdef ONLXPLUS
CC = /cvmfs/sft.cern.ch/lcg/releases/gcc/9.2.0/x86_64-centos7/bin/g++
else
CC = g++
endif

CFLAGS = -Wall -O3 -std=c++11 -Wno-narrowing
LIBS   = -lz -lm -lstdc++fs -pthread

PROG = ../bin/analysis

OBJS = ProcessData.o Physics.o Distributions.o Selection.o Luminosity.o \
       TrkReco.o TrkEffic.o Resolutions.o ParticleId.o MostProbable.o \
       RpReco.o RpSimu.o RpVeto.o RpEffic.o RpFit.o RpPat.o RpPrint.o RpPolygon.o \
       SimplexMinimizer.o HistoVal.o Histo.o \
       RandomSobol.o Random.o Helper.o gzstream.o

all : $(PROG)

#
../bin/analysis : analysis.o $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

objs : $(OBJS)
	@ echo all objects compiled

%.o : %.cc $(wildcard ../interface/*.h)
	$(CC) $(CFLAGS) -c -o $@ $<

#
clean :
	rm -f *.o $(PROG)
