# Include files
SOURCES=initLB.c visualLB.c boundary.c collision.c streaming.c computeCellValues.c main.c helper.c

# Compiler
# --------
CC=gcc

CFLAGS=-Werror -pedantic -Wall -std=c99 -g #-enable-checking -O0 -v -da -Q

# Linker flags
# ------------
LDFLAGS= -lm -pg

OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE=lbsim

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS)

clean:
	rm -f $(OBJECTS) $(EXECUTABLE)

clean_debug:
	rm -f *.expand *.vregs *.into_cfglayout *.jump *.reginfo *.outof_cfglayout *.split1 *.dfinit *.mode_sw *.asmcons *.ira *.reload *.split2 *.pro_and_epilogue *.jump2 *.stack *.alignments *.mach *.barriers *.shorten *.nothrow *.dwarf2 *.final *.dfinish


$(OBJECTS): %.o : %.c
	$(CC) $(CFLAGS) -c $< -o $@
