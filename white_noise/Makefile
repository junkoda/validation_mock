#
# Isotropic Distortion of Gaussian field
#

EXEC      = white_mock
DIRS      = $(GSL_DIR)

all: $(EXEC)

INCLDIRS  = $(foreach dir, $(DIRS), -I$(dir)/include)
LIBDIRS   = $(foreach dir, $(DIRS), -L$(dir)/lib) 
CXXFLAGS := $(INCLDIRS)

OBJS := white_mock.o

LIBS := $(LIBDIRS) -lm -lboost_program_options
LIBS += -lgsl -lgslcblas

white_mock: $(OBJS)
	$(CXX) $(OBJS) $(LIBS) -o $@


.PHONY: clean
clean :
	rm -f $(EXEC) $(OBJS)

