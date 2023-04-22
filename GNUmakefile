##################################
## Elliptica ID reader Makefile ##
##################################

## lib name
LIBNAME = elliptica_id_reader

## compiler settings ('?=' means if not set already)
CC ?= gcc
CFLAGS ?= -fopenmp -O3 -std=c99 -pedantic
LIB_CFLAGS = -fPIC 

## Creating a static or dynamic library. Note you can make both lib. types by 
## calling make twice; once with the default LIB_TYPE and then with the other option,
## for instance, $ make LIB_TYPE="dynamic".
## This variable can be useful when compiling in MacOS where we want to use static 
## and not ".dylib" lib.
LIB_TYPE ?= "static" ## "dynamic"

## ---------------------------------------------------------------------- ##
## ---------------------------------------------------------------------- ##

BASE=$(shell /bin/pwd)
SRCD=$(BASE)/src
INCD=$(BASE)/include
OBJD=$(BASE)/obj
LIBD=$(BASE)/lib

SRC=$(wildcard $(SRCD)/*.c)
OBJ=$(patsubst $(SRCD)/%.c,$(OBJD)/%.o,$(SRC))
INC_PARAMS=-I$(SRCD)

## a quick setup of lib
ifeq (dynamic,$(strip $(LIB_TYPE)))
	LIB=$(LIBD)/lib$(LIBNAME).so
else
	LIB=$(LIBD)/lib$(LIBNAME).a
endif


## archive cmd
AR = ar

.PHONY: all clean

all: $(LIB)
	@echo '---'
	@echo "'$(LIBNAME)' library made successfully :)"
	@echo '---'

$(OBJD)/%.o: $(SRCD)/%.c
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(LIB_CFLAGS) $(INC_PARAMS) -c $< -o $@

## only a dynamic lib ##
$(LIBD)/lib$(LIBNAME).so: $(OBJ)
	@echo '---'
	@echo "Making a dynamic library ... "
	@echo '---'
## cp header
	@mkdir -p $(INCD)
	@cp $(SRCD)/elliptica_id_reader_lib.h $(INCD)
## shared lib
	@mkdir -p $(LIBD)
	$(CC) $(CFLAGS) $(LIB_CFLAGS) -shared $(OBJ) -o $(LIB)

## only a static lib ##
$(LIBD)/lib$(LIBNAME).a: $(OBJ)
	@echo '---'
	@echo "Making a static library ... "
	@echo '---'
## cp header
	@mkdir -p $(INCD)
	@cp $(SRCD)/elliptica_id_reader_lib.h $(INCD)
## static lib
	@mkdir -p $(LIBD)
	$(AR) rcs $(LIB) $(OBJ)

## remove installation
clean:
	@echo '---'
	@echo "Removing ..."
	@rm -rf $(OBJD)
	@rm -rf $(LIBD)
	@rm -rf $(INCD)
	@echo "Done!"
	@echo '---'
