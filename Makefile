CC=gcc
LD=${CC} -Xlinker -Map=.program.map -Xlinker --cref 
CFLAGS= -Wall -pedantic -std=c99 -O3 -g -DLOCAL -DSORTEDUNMAPPED -D_LARGEFILE_SOURCE  -DFILEBUFFEREDMERGE -D_FILE_OFFSET_BITS=64 -DDBGNFO -DSHOWALIGN -DDBGLEVEL=0 -DPROGNFO -Ilibs -Ilibs/sufarray
INC := -I include
CTAGS = ctags > tags
LIB = -lm -lpthread -lz -lncurses -L/usr/local/lib -L libs -lform -lmenu -lgmp -lgsl -lgslcblas


LIBDIR := libs
BUILDDIR:= build
TARGETDIR := .
TARGETEXT := .x
SRCEXT := c
PRGTARGETS := hydi filebuffer_test

SOURCES := $(shell find $(LIBDIR) -type f -name *.$(SRCEXT))
EXCLUDE := $(LIBDIR)/splines.c $(LIBDIR)/snvsplines.c
SOURCES := $(filter-out $(EXCLUDE), $(SOURCES)) 

PRGSOURCES := $(patsubst %,$(LIBDIR)/%.c,$(PRGTARGETS))
LIBSOURCES := $(filter-out $(PRGSOURCES), $(SOURCES)) 

OBJECTS := $(patsubst $(LIBDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
LIBOBJECTS :=  $(patsubst $(LIBDIR)/%,$(BUILDDIR)/%,$(LIBSOURCES:.$(SRCEXT)=.o))


$(PRGTARGETS): $(OBJECTS)
	@echo "Linking $@";	
	$(CC) $(LIBOBJECTS) $(BUILDDIR)/$@.o -o $(TARGETDIR)/$@$(TARGETEXT) $(LIB)


$(BUILDDIR)/%.o: $(LIBDIR)/%.c
	@echo "Building... library for source $@";
	@mkdir -p $(BUILDDIR)	
	$(CC) $(CFLAGS) $(INC) -c -o $@ $<

clean:
	@echo " Cleaning...";
	@echo " $(RM) -r $(BUILDDIR) $(PRGTARGETS)"; $(RM) -r $(BUILDDIR) $(PRGTARGETS)

all: $(PRGTARGETS)
.PHONY: clean
