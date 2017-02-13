#
# Isotropic Distortion of Gaussian field
#

EXEC      = white_mock white_power
DIRS      = $(GSL_DIR)

all: $(EXEC)

INCLDIRS  = $(foreach dir, $(DIRS), -I$(dir)/include)
LIBDIRS   = $(foreach dir, $(DIRS), -L$(dir)/lib) 
CXXFLAGS := $(INCLDIRS)

OBJS1 := white_mock.o

LIBS1 := $(LIBDIRS) -lm -lboost_program_options
LIBS1 += -lgsl -lgslcblas

white_mock: $(OBJS1)
	$(CXX) $(OBJS1) $(LIBS) -o $@


OBJS2 := white_power.o
LIBS2 := $(LIBDIRS) -lm -lboost_program_options

white_power: $(OBJS2)
	$(CXX) $(OBJS2) $(LIBS2) -o $@




.PHONY: clean
clean :
	rm -f $(EXEC) $(OBJS)

