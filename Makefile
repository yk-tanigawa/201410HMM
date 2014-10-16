CC = g++
LD = g++
CFLAGS =  -O2
LDFLAGS =  #-lpthread
SRCS := $(wildcard *.cpp) # wildcard
OBJS = $(SRCS:.cpp=.o)
DEPS = $(SRCS:.cpp=.dep)
EXEC = $(SRCS:.cpp=)
RM = rm -f


all:

kadai1main: kadai1main.o
	$(LD)  -o $@ $^ $(LDFLAGS)

kadai2main: kadai2main.o
	$(LD)  -o $@ $^ $(LDFLAGS)

kadai3main: kadai3main.o
	$(LD)  -o $@ $^ $(LDFLAGS)

clean:
	$(RM) $(OBJS) $(EXEC) *~

check-syntax:
	$(CC) -Wall -Wextra -pedantic -fsyntax-only $(CHK_SOURCES)

.PHONY: clean 
	all clean

