# copyright: Michael Safyan

program_NAME := cc.bin
program_C_SRCS := $(wildcard *.c)
program_CXX_SRCS := $(wildcard *.cpp)
program_C_OBJS := ${program_C_SRCS:.c=.o}
program_CXX_OBJS := ${program_CXX_SRCS:.cpp=.o}
program_OBJS := $(program_C_OBJS) $(program_CXX_OBJS)
program_INCLUDE_DIRS := ${BOOST_ROOT} 
program_LIBRARY_DIRS := 
program_LIBRARIES :=

CPPFLAGS += -std=c++0x  -O3 -D__USE_XOPEN2K8 -fopenmp -DSINGLEPRECISION  
CPPFLAGS += $(foreach includedir,$(program_INCLUDE_DIRS),-I$(includedir))
LDFLAGS +=  -pthread 
LDFLAGS += $(foreach library,$(program_LIBRARIES),-l$(library))
LDFLAGS += -lz
.PHONY: all clean distclean

all: $(program_NAME)

$(program_NAME): $(program_OBJS)
	$(LINK.cc) $(program_OBJS) -o $(program_NAME) $(LDFLAGS)

test: test/log/check.txt
	./$(program_NAME) -v

test/log/check.txt: $(program_NAME)
	mkdir -p test/log/
	@bash test.sh 

clean: cleantest
	@- $(RM) $(program_NAME)
	@- $(RM) $(program_OBJS)

cleantest:
	@- $(RM) test/out/*
	@- $(RM) test/out2/*
	@- $(RM) test/log/check.txt

distclean: clean

