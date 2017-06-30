# =======================================
#   Makefile for MpdRoot Analysis team
#     by V. Kireyeu and A. Mudrokh
#             v2017-04-18
#
#  This Makefile uses for analysis script 
#  compilation via the system compiler.
#  Please do not forget to run config.sh
#  script from MpdRoot build directory!!! 
# =======================================

# Used compiler
CC=g++

# Makefile will proceed all files with .cc suffix
# SOURCES=$(wildcard src/*.cc)
# Output executable has the same name without suffix
# OBJECTS=$(patsubst src/%.cc,%, $(SOURCES))

SOURCES=$(wildcard *.cc)
OBJECTS=$(patsubst %.cc,%, $(SOURCES))

ROOTCONFIG := root-config
# Path to the MpdRoot installation
MPDROOTDIR := /home/alex/git/march_2017/mpdroot


CFLAGS := $(shell $(ROOTCONFIG) --cflags)
CFLAGS += -I$(MPDROOTDIR)/base/event 
CFLAGS += -I$(MPDROOTDIR)/base/source 
CFLAGS += -I$(MPDROOTDIR)/base/steer 
CFLAGS += -I$(MPDROOTDIR)/base/field 
CFLAGS += -I$(MPDROOTDIR)/base/sim 
CFLAGS += -I$(MPDROOTDIR)/geobase 
CFLAGS += -I$(MPDROOTDIR)/tpc 
CFLAGS += -I$(MPDROOTDIR)/kalman 
CFLAGS += -I$(MPDROOTDIR)/lhetrack 
CFLAGS += -I$(MPDROOTDIR)/mcstack 
CFLAGS += -I$(MPDROOTDIR)/strawendcap 
CFLAGS += -I$(MPDROOTDIR)/etof 
CFLAGS += -I$(MPDROOTDIR)/tof 
CFLAGS += -I$(MPDROOTDIR)/sft 
CFLAGS += -I$(MPDROOTDIR)/sts 
CFLAGS += -I$(MPDROOTDIR)/parbase 
CFLAGS += -I$(MPDROOTDIR)/mpddata 
CFLAGS += -I$(MPDROOTDIR)/mpdbase 
CFLAGS += -I$(MPDROOTDIR)/fairtools 
CFLAGS += -I$(MPDROOTDIR)/clustering 
CFLAGS += -I$(MPDROOTDIR)/geometry 
CFLAGS += -I$(MPDROOTDIR)/mpdpid 

LDFLAGS := $(shell $(ROOTCONFIG) --ldflags)
LDFLAGS += $(shell $(ROOTCONFIG) --glibs)
LDFLAGS += -L$(MPDROOTDIR)/build/lib -lMpdField -lFairTools -lGeoBase 
LDFLAGS += -lBase -lMCStack -lMpdField -lPassive -lTrkBase 
LDFLAGS += -lMpdBase -lMpdData -lMpdPid -lMpdGeneralGenerator 
LDFLAGS += -lKalman -lGeane -lCluster -ltpc -lTof -lEtof -lStrawendcap 
LDFLAGS += -lZdc -lFfd -lCpc -lStrawECT -lLHETrack -lSts -lCluster 

all:	$(OBJECTS)

# $(OBJECTS): % : src/%.cc
$(OBJECTS): % : %.cc
	$(CC) $(CFLAGS) $(LDFLAGS) $< -o $@

clean:
	rm -vf $(OBJECTS)