DIRS = lib exec

all: $(DIRS)

.PHONY: lib exec clean

#
# Set compilers and options
#
CC      := mpic++ -std=c++11
CXX     := mpic++ -std=c++11
OPT     := -DDOUBLEPRECISION

export CC CXX OPT


lib:
	cd lib && $(MAKE) lib

exec:
	cd exec && $(MAKE) exec

clean:
	for dir in $(DIRS); do (cd $$dir && $(MAKE) clean); done
