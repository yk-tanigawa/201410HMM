CC = g++
LD = g++
CFLAGS = -Wall -lm  #-O2
LDFLAGS = -lm #-lpthread
SRCS := $(wildcard *.cpp) # wildcard
OBJS = $(SRCS:.cpp=.o)
DEPS = $(SRCS:.cpp=.dep)
EXEC = $(SRCS:.cpp=)
RM = rm -f


all: hmm

hmm: hmm.o
	$(LD)  -o $@ $^ $(LDFLAGS)

clean:
	$(RM) $(OBJS) $(EXEC) *~

.PHONY:
	all clean
