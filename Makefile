#    Copyright (C) 2011-2104 University of Southern California and
#                            Andrew D. Smith and Timothy Daley
#
#    Authors: Timothy Daley and Andrew D. Smith
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
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#


ifndef ROOT
ROOT = $(shell pwd)
endif

ifndef SMITHLAB_CPP
SMITHLAB_CPP=$(ROOT)/smithlab_cpp/
endif

SOURCES = $(wildcard *.cpp)
OBJECTS = $(patsubst %.cpp,%.o,$(SOURCES))
PROGS = SamplePoissonDirichlet
INCLUDEDIRS = $(SMITHLAB_CPP) 
INCLUDEARGS = $(addprefix -I,$(INCLUDEDIRS))

LIBS += -lgsl -lgslcblas -lz

CXX = g++ 
CXXFLAGS = -Wall -fPIC -fmessage-length=50

# Flags for Mavericks
ifeq "$(shell uname)" "Darwin"
CXXFLAGS += -arch x86_64
ifeq "$(shell if [ `sysctl -n kern.osrelease | cut -d . -f 1` -ge 13 ];\
              then echo 'true'; fi)" "true"
CXXFLAGS += -stdlib=libstdc++
endif
endif


OPTFLAGS = -O2
DEBUGFLAGS = -g -lefence -lpthread -L/usr/local/lib/

ifdef DEBUG
CXXFLAGS += $(DEBUGFLAGS)
endif



ifdef OPT
CXXFLAGS += $(OPTFLAGS)
endif

all: $(PROGS)

$(PROGS): $(addprefix $(SMITHLAB_CPP)/, \
          smithlab_os.o smithlab_utils.o OptionParser.o)

SamplePoissonDirichlet: PoissonDirichlet.o

%.o: %.cpp %.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(INCLUDEARGS)

%: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^ $(INCLUDEARGS) $(LIBS)

install: $(PROGS)
	@mkdir -p $(ROOT)/bin
	@install -m 755 $(PROGS) $(ROOT)/bin

clean:
	@-rm -f $(PROGS) *.o *~

.PHONY: clean
