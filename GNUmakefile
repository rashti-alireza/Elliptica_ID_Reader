# Makefile for Elliptica

BASE=$(shell /bin/pwd)
SRCD=$(BASE)/src
##INCD=$(BASE)/include
INCD=$(BASE)/src
OBJD=$(BASE)/obj
LIBD=$(BASE)/lib


NAME=EllipticaRun
LIBNAME=elliptica_id_reader
EXE=$(NAME).x
SRC=$(wildcard $(SRCD)/*.c)
OBJ=$(patsubst $(SRCD)/%.c,$(OBJD)/%.o,$(SRC))
INC=$(wildcard $(INCD)/*.c)
INC_PARAMS=$(foreach d, $(INCD), -I$d)

OBJ2=$(wildcard $(OBJD)/*.o)
LIB=$(LIBD)/lib$(LIBNAME).so


# mandatory flags

CC = gcc
LD = ld
AR = ar

CFLAGS = -std=c99 -fPIC -pedantic
CFLAGS+= -I/usr/include/suitesparse
LDLFLAGS=-lumfpack -lm

#CFLAGS += -O3
CFLAGS += -g

all: $(LIB)
	@echo "All done"

$(EXE): $(OBJ)
	@echo "Building $@ ..."
	@echo

	$(CC) $(CFLAGS) $(LFLAGS) $(OBJ) -o $@ $(LDLFLAGS)
        
#	# this doesn't need to be a shared object
#	@rm $(OBJD)/EllipticaRun.o > /dev/null 2>&1

$(OBJD)/%.o: $(SRCD)/%.c
	@echo "Building objects ..."
	@echo

	@mkdir -p $(dir $@)
#	@echo "Compiling $< ... -> $@"
	$(CC) $(CFLAGS) $(INC_PARAMS) -c $< -o $@

$(LIB): $(OBJ)
	@echo "Making libraries... "
	@echo

	@mkdir -p $(LIBD)
	$(CC) $(CFLAGS) -shared $(OBJ) -o $(LIBD)/lib$(LIBNAME).so

	$(AR) rcs $(LIBD)/libElliptica.a $(OBJ)
#	$(AR) rcs $(LIBD)/libElliptica_static.a $@ $^

clean:
	@echo "Cleaning ..."
	@rm -rf $(OBJD)
	@rm -rf $(LIBD)
	@rm -rf $(EXE)
	@echo "... done"

.PHONY: all clean
