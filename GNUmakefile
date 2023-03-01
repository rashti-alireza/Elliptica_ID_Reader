##################################
## Elliptica ID reader Makefile ##
##################################

## lib name
LIBNAME = elliptica_id_reader

## compiler settings
CC = gcc

CFLAGS = -std=c99 -fPIC -pedantic
CFLAGS += -fopenmp
#CFLAGS += -O3

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
LIB=$(LIBD)/lib$(LIBNAME).so

## archive cmd
AR = ar

.PHONY: all clean

all: $(LIB)
	@echo '---'
	@echo "'$(LIBNAME)' library made successfully :)"
	@echo '---'

$(OBJD)/%.o: $(SRCD)/%.c
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(INC_PARAMS) -c $< -o $@

$(LIB): $(OBJ)
	@echo '---'
	@echo "Making libraries ... "
	@echo '---'

	@mkdir -p $(LIBD)
	
## shared lib
	$(CC) $(CFLAGS) -shared $(OBJ) -o $(LIBD)/lib$(LIBNAME).so

## static lib
	$(AR) rcs $(LIBD)/lib$(LIBNAME).a $(OBJ)

## cp header
	@mkdir -p $(INCD)
	@cp -u $(SRCD)/elliptica_id_reader_lib.h $(INCD)
	
clean:
	@echo '---'
	@echo "Removing ..."
	@rm -rf $(OBJD)
	@rm -rf $(LIBD)
	@rm -rf $(INCD)
	@echo "Done!"
	@echo '---'
