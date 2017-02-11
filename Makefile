DIRS = lib

.PHONY: lib clean check
all: $(DIRS)

#
# Set compilers and options
#
CC      := mpic++ -std=c++11
CXX     := mpic++ -std=c++11
OPT     := -DDOUBLEPRECISION

export CC CXX OPT


lib:
	cd lib && $(MAKE) lib

clean:
	for dir in $(DIRS); do (cd $$dir && $(MAKE) clean); done
