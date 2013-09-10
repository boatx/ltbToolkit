SHELL = /bin/sh
.SUFFIXES:
.SUFFIXES: .c .o .d .a
src = ../../src/
lib = ../../src/.libs/libltb.a
#lib = /home/lodka/temp/lib/usr/local/lib/ltbtoolkit/libltb.so
VPATH = src:../../src
BINDIR = bin
OBJS = flrw.o main.o
CC = gcc
DEBUG = -g
CFLAGS = -Wall -Wextra -pedantic -std=c99 -O1 -fopenmp
LDFLAGS = -lgsl -lgslcblas -lm -lgomp

$(BINDIR)/flrw.out: $(OBJS) $(lib)  | $(BINDIR)
	$(CC) $^ $(DEBUG) -o $@ $(LDFLAGS) 

all : $(BINDIR)/ltb.out
.PHONY: all

%.d : $(OBJS:.o=.c)
	@set -e; rm -f $@; \
	$(CC) -MM $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$
-include $(OBJS:.o=.d)

$(OBJS): %.o: %.c
	$(CC) -c -I$(src) $(CFLAGS) $(DEBUG) $< -o $@

clean:
	rm -rf *.o *.d  bin/
.PHONY : clean

cleandat:
	rm -rf *.dat
.PHONY : cleandat

$(BINDIR) :
	mkdir $(BINDIR)

