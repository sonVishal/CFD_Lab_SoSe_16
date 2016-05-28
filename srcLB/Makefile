# Include files
SOURCES=initLB.c visualLB.c boundary.c collision.c streaming.c computeCellValues.c main.c helper.c

# Compiler
# --------
CC=gcc

CFLAGS=-Werror -g -pedantic -Wall -std=c99 #-O0 -g -enable-checking -v -da -Q

# Linker flags
# ------------
LDFLAGS= -lm #-pg

OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE=lbsim

all: CFLAGS+= -O0 
all: clean $(EXECUTABLE) ; @echo "\nINFO:Use make speed for a faster code without checks."

speed: CFLAGS+= -DNDEBUG -O3
speed: clean $(EXECUTABLE)

debug: CFLAGS+=-O0 -g -enable-checking -v -da -Q
debug: LDFLAGS+=-pg
debug: clean $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS)

clean:
	rm -f $(OBJECTS) $(EXECUTABLE)

clean_debug:
	rm -f *.expand *.vregs *.into_cfglayout *.jump *.reginfo *.outof_cfglayout *.split1 *.dfinit *.mode_sw *.asmcons *.ira *.reload *.split2 *.pro_and_epilogue *.jump2 *.stack *.alignments *.mach *.barriers *.shorten *.nothrow *.dwarf2 *.final *.dfinish


$(OBJECTS): %.o : %.c
	$(CC) $(CFLAGS) -c $< -o $@
