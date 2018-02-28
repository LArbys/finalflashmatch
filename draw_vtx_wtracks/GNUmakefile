
# include options for this package
INCFLAGS  = -I.
INCFLAGS += $(shell root-config --cflags)
INCFLAGS += $(shell larlite-config --includes)
INCFLAGS += -I$(LARLITE_USERDEVDIR)
INCFLAGS += $(shell larcv-config --includes)
INCFLAGS += $(shell larlitecv-config --includes)

LDFLAGS += $(shell root-config --ldflags --libs)
LDFLAGS += $(shell larlite-config --libs)
LDFLAGS += $(shell larcv-config --libs)
LDFLAGS += $(shell larlitecv-config --libs)
#LDFLAGS += -lLArOpenCV_Core
#LDFLAGS += -lLArOpenCV_ImageClusterBase
#LDFLAGS += -lLArOpenCV_ImageClusterAlgoFunction
#LDFLAGS += -lLArOpenCV_ImageClusterAlgoData
#LDFLAGS += -lLArOpenCV_ImageClusterAlgoClass

# platform-specific options
OSNAME          = $(shell uname -s)
HOST            = $(shell uname -n)
OSNAMEMODE      = $(OSNAME)



all: draw_vtx_wtracks

draw_vtx_wtracks: draw_vtx_wtracks.cxx
	g++ -g -fPIC $(INCFLAGS) draw_vtx_wtracks.cxx -o draw_vtx_wtracks $(LDFLAGS)

