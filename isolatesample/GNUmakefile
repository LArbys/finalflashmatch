
# include options for this package
INCFLAGS  = -I.
INCFLAGS += $(shell root-config --cflags)
INCFLAGS += $(shell larlite-config --includes)
INCFLAGS += -I$(LARLITE_USERDEVDIR)
INCFLAGS += $(shell larcv-config --includes)
INCFLAGS += $(shell larlitecv-config --includes)
INCFLAGS += -I$(LAROPENCV_BASEDIR)
INCFLAGS += -I$(LARCV_INCDIR)/LArOpenCVHandle/

LDFLAGS += $(shell root-config --ldflags --libs)
LDFLAGS += $(shell larlite-config --libs)
LDFLAGS += $(shell larcv-config --libs)
LDFLAGS += $(shell larlitecv-config --libs)

LDFLAGS += -lLArOpenCV_Core
LDFLAGS += -lLArOpenCV_ImageClusterBase
LDFLAGS += -lLArOpenCV_ImageClusterAlgoFunction
LDFLAGS += -lLArOpenCV_ImageClusterAlgoData
LDFLAGS += -lLArOpenCV_ImageClusterAlgoClass

# note: llcvprocessor headers and libraries are in larlitecv/build/include and lib
#LDFLAGS += -l

# platform-specific options
OSNAME          = $(shell uname -s)
HOST            = $(shell uname -n)
OSNAMEMODE      = $(OSNAME)



all: select_tune_sample

select_tune_sample: select_tune_sample.cxx
	g++ -g -fPIC $(INCFLAGS) select_tune_sample.cxx -o select_tune_sample $(LDFLAGS)

