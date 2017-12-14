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
#DONT FUCK AROUND WITH SRCDIR, IT DOESN'T WORK

NDIM=3

# Define compile time variables to be passed on
DEFINES=-DNDIM=$(NDIM)

#==================================
# COMPILER OPTIONS AND DIRECTIVES
#==================================


#Compiler to use
CC=gcc -I$(IDIR) $(DEFINES)

# flags for gfortran compiler
CFLAGS= -g -O0 -Wall # debug
#CFLAGS= -O3  		 # performance



#==================================
# FILE LISTS
#==================================

#include paths. Will be followed in that order.
VPATH = $(IDIR):../src

#SOURCES=$(wildcard *.c)

SOURCES= commons.c  init.c  main_program.c

HEADERS = commons.h  init.h

OBJECTS = $(notdir $(patsubst %.c, %.o, $(wildcard ../src/*.c)))



#==================================
# RECEPIES
#==================================


%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c -o $(notdir $@) $<
	$(info $@)
#	-c : compile and assemble, but do not link

$(EXEC): $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $^



#==================================
#==================================

.PHONY: clean # keeps makefile from doing something with a file named clean

clean:
	rm -f *.o $(EXEC)

