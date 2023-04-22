## a quick makefile for test purposes

## set paths for each compiling option of interest
INC_DIRS_V1 += /home/alireza/Workstation/Elliptica/ID_Reader/include
LIB_DIRS_V1 += /home/alireza/Workstation/Elliptica/ID_Reader/lib

INC_DIRS_V2 += /home/alireza/Workstation/Elliptica_ID_Reader/include
LIB_DIRS_V2 += /home/alireza/Workstation/Elliptica_ID_Reader/lib

SRC = id_reader_nsns.c
EXE = id_reader_nsns
CC ?= gcc
CFLAGS ?= -fopenmp -O3

LIB_NAME = -lelliptica_id_reader

.PHONY: clean from_code from_reader

## use id reader lib from the elliptica code we further need umfpack
from_code:
	$(CC) -o $(EXE) $(CFLAGS) $(SRC) -I$(INC_DIRS_V1) -L$(LIB_DIRS_V1) \
	$(LIB_NAME) -Wl,-rpath,$(LIB_DIRS_V1) -lumfpack -lm
	
 

## use id reader lib from the stand alone Elliptica_ID_Reader repo
from_reader:
	$(CC) -o $(EXE) $(CFLAGS) $(SRC) -I$(INC_DIRS_V2) -L$(LIB_DIRS_V2) \
	$(LIB_NAME) -Wl,-rpath,$(LIB_DIRS_V2) -lm
	

clean:
	@rm -rf $(EXE)
 
