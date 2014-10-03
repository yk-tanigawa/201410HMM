CC = gcc
LD = gcc
CFLAGS = -Wall -lm  #-O2
LDFLAGS = -lm #-lpthread
SRCS := $(wildcard *.c) # wildcard
OBJS = $(SRCS:.c=.o)
DEPS = $(SRCS:.c=.dep)
EXEC = $(SRCS:.c=)
RM = rm -f


all: viterbi

viterbi.o: viterbi.c viterbi.h

viterbi: viterbi.o hmm.o
	$(LD)  -o $@ $^ $(LDFLAGS)

clean:
	$(RM) $(OBJS) $(EXEC) *~

.PHONY:
	all clean
