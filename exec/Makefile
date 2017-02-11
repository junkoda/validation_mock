EXEC       = mock_gaussian

default:
	cd .. && $(MAKE) exec

exec: $(EXEC)

# Define OPENMP to enable MPI+OpenMP hybrid parallelization
# OPENMP  = -fopenmp # -openmp for Intel, -fopenmp for gcc, llvm doesn't support

# Library location if not in a standard path
FFTW3_DIR ?= #e.g. /Users/jkoda/Research/opt/gcc/fftw3
GSL_DIR   ?= #e.g. /Users/jkoda/Research/opt/gcc/gsl

CPPFLAGS  := -O3 $(OPENMP) $(OPT) -I../lib


DIR_PATH   = $(FFTW3_DIR) $(GSL_DIR)

CPPFLAGS  += $(foreach dir, $(DIR_PATH), -I$(dir)/include)
LIBS      += $(foreach dir, $(DIR_PATH), -L$(dir)/lib)

LIBS      += -lm -L../lib -lfs
LIBS      += -lgsl -lgslcblas -lboost_program_options

ifeq (,$(findstring -DDOUBLEPRECISION, $(OPT)))
  # Single precision FFTW
  FFTWSUF=f
endif
LIBS += -lfftw3$(FFTWSUF) -lfftw3$(FFTWSUF)_mpi

ifdef OPENMP
  LIBS += -lfftw3$(FFTWSUF)_omp
  #LIBS += -lfftw3$(FFTWSUF)_threads # for thread parallelization instead of omp
endif

OBJS = mock_gaussian.o

mock_gaussian: $(OBJS)
	$(CXX) $(OBJS) $(LIBS) -o $@

.PHONY: clean
clean :
	rm -f $(EXEC) $(OBJS)
