CC = g++
CFLAGS	= -I$(HOME)/usr/include
LDFLAGS	= -L$(HOME)/usr/lib
LDLIBS	= -lfftw3_threads -lfftw3 -lm -fopenmp

SRCDIR = src
RUNDIR = run
DATADIR = data
EXECUTABLE = AAA 

SOURCES	= $(wildcard $(SRCDIR)/*.cpp)
OBJECTS	= $(patsubst $(SRCDIR)/%.cpp,$(RUNDIR)/%.o,$(SOURCES))

export LD_RUN_PATH=$(HOME)/usr/lib

all: $(RUNDIR)/$(EXECUTABLE)

$(RUNDIR)/$(EXECUTABLE): $(OBJECTS)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^ ${LDLIBS}

$(OBJECTS): $(RUNDIR)/%.o : $(SRCDIR)/%.cpp
	$(CC) $(CFLAGS) $(LDFLAGS) -c $< ${LDLIBS} -o $@ 

clean:
	rm -rf $(RUNDIR)/*.o
	rm -rf $(RUNDIR)/$(EXECUTABLE)
	rm -rf qsub.sh.* 
	rm -rf $(DATADIR)/*.vtk
	rm -rf $(DATADIR)/*.csv

