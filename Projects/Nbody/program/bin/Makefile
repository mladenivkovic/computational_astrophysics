# This Makefile is used to compile the main_program.c code correctly.
#
# SYNTAX:
#target: depencencies
#[tab] 	system command

#==================================
# COMPILATION TIME PARAMETERS
#==================================

EXEC=nbody #executable name
IDIR=../include #include headers from this directory
SRCDIR=../src


#==================================
# COMPILER OPTIONS AND DIRECTIVES
#==================================


#Compiler to use
CC=gcc -I$(IDIR) -fopenmp



# flags for gcc compiler
CFLAGS= -g -Og -Wall -Wextra     #debug
# CFLAGS= -O3  		 # performance


LDFLAGS= -lm


#==================================
# FILE LISTS
#==================================

#include paths. Will be followed in that order.

VPATH=$(IDIR):$(SRCDIR)


OBJECTS = $(notdir $(patsubst %.c, %.o, $(wildcard $(SRCDIR)/*.c)))



#==================================
# RECEPIES
#==================================


default: $(EXEC)


$(EXEC): $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

%.o:  %.c %.h
	$(CC) $(CFLAGS) -o $@ -c $<
#-c : compile and assemble, but do not link



#==================================
#==================================

.PHONY: clean # keeps makefile from doing something with a file named clean

clean:
	rm -f $(OBJECTS) $(EXEC)

