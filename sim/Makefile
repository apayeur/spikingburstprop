# shorthand for your MPI C++ compiler
CC = mpicxx
# Your Auryn install path
AURYNDIR = $(HOME)/auryn
# Path to libauryn.a (see compiling Auryn)
BUILDDIR = $(AURYNDIR)/build/release/

CFLAGS= -ansi -Wall -pipe -O3 -DNDEBUG -std=c++11 -ffast-math -funsafe-math-optimizations \
        -march=native -mtune=native -pedantic \
        -I/opt/local/include -I../src -I$(AURYNDIR)/src/auryn -I$(AURYNDIR)/src

LDFLAGS=$(BUILDDIR)/src/libauryn.a \
		-L/opt/local/lib -lboost_program_options-mt -lboost_serialization-mt -lboost_mpi-mt -lboost_system-mt -lboost_filesystem-mt

# debug
# BUILDDIR = $(AURYNDIR)/build/debug/
# CFLAGS= -ansi -Wall -pipe -g \
#         -march=native -mtune=native -pedantic \
#         -I$(AURYNDIR)/src -I$(AURYNDIR)/dev/src


OBJECTS=STPeTMConnection.o BCPConnection.o EBCPConnection.o SineCurrentInjector.o 
#TransmitBurstConnection.o TransmitEventConnection.o 

SIMULATIONS=sim_plasticity_rule sim_propagation sim_propagation_ficurve sim_noise_matching

.SECONDARY:

all: $(SIMULATIONS) 

sim_%: sim_%.o $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) $< $(LDFLAGS) -o $(subst .o,,$<)

%.o : %.cpp
	$(CC) $(CFLAGS) -c $<

%.o : ../src/%.cpp
	$(CC) $(CFLAGS) -c $<

#figs: 

#plot_%.tex: plot_%.gnu 
#	gnuplot $< 
#	epstopdf epslatex-inc.eps
#	mv epslatex.out $@

#%.pdf: %.tex
#	pdflatex $< 

clean: 
	rm -f *.o $(SIMULATIONS)
