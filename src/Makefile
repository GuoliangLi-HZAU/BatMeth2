#    Copyright (C) 2013 University of Southern California and
#                       Egor Dolzhenko
#                       Andrew D Smith
#
#    Authors: Andrew D. Smith and Egor Dolzhenko
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#ifndef CPP
#$(error CPP variable undefined)
#endif

CPP = $(shell pwd)
SOURCES = $(wildcard *.cpp)
OBJECTS = $(patsubst %.cpp,%.o,$(SOURCES))
PROGS = MOMO

CXX = g++
CFLAGS = -Wall -fmessage-length=50
CPPFLAGS =
CXXFLAGS = -Wall -fmessage-length=50
OPTFLAGS = -O2
DEBUGFLAGS = -g
CXXFLAGS =

ifdef DEBUG
CXXFLAGS += $(DEBUGFLAGS)
endif

ifdef OPT
CXXFLAGS += $(OPTFLAGS)
endif

COMMON_DIR = ./ # ../common
INCLUDEDIRS =  $(CPP) $(COMMON_DIR)
INCLUDEARGS = $(addprefix -I,$(INCLUDEDIRS))

LIBS = -lz

# Flags for Mavericks
ifeq "$(shell uname)" "Darwin"
CXXFLAGS += -arch x86_64
INCLUDEARGS += -I/usr/local/include
LIBS += -L/usr/local/lib
ifeq "$(shell if [ `sysctl -n kern.osrelease | cut -d . -f 1` -ge 13 ];\
							then echo 'true'; fi)" "true"
CXXFLAGS += -stdlib=libstdc++
endif
endif

nogzip: $(PROGS)

all: $(PROGS) 

install: $(PROGS)
	@mkdir -p $(SRC_ROOT)/bin
	@install -m 755 $(PROGS) $(SRC_ROOT)/bin

OBJS= regression.o combine_pvals.o merge.o

MOMO:
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(CXXFLAGS) -o genome2cg genome2cg.cpp $(INCLUDEARGS) $(LIBS)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(CXXFLAGS) -o genomebinLen genomebinLen.cpp $(INCLUDEARGS) $(LIBS)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(CXXFLAGS) -o memalign memalign.cpp $(INCLUDEARGS) -I./samtools-0.1.18/ -L./samtools-0.1.18/ -lbam -lz -lpthread $(LIBS)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(CXXFLAGS) -o build_index_rrbs build_index_rrbs.cpp $(INCLUDEARGS) $(LIBS)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(CXXFLAGS) calmeth.cpp -o calmeth -m64 -I./samtools-0.1.18/ -L./samtools-0.1.18/ -lbam -lz
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(CXXFLAGS) splitSam.cpp -o splitSam -m64 -I./samtools-0.1.18/ -L./samtools-0.1.18/ -lbam -lz -pthread $(INCLUDEARGS) $(LIBS)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(CXXFLAGS) -o methyGff methyGff.cpp $(INCLUDEARGS) $(LIBS)

%.o : %.cpp %.hpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $< $(INCLUDEARGS) $(LIBS)

%: %.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(CXXFLAGS) -o $@ $^ $(INCLUDEARGS) $(LIBS)

clean:
	@-rm -f $(PROGS) *.o *.so *.a *~ genome2cg genomebinLen memalign build_index_rrbs calmeth splitSam
