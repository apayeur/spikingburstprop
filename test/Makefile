# shorthand for your MPI C++ compiler
CC = mpicxx
# Your Auryn install path
AURYNDIR = $(HOME)/auryn
# Path to libauryn.a (see compiling Auryn)
BUILDDIR = $(AURYNDIR)/build/release/

CFLAGS= -ansi -Wall -pipe -O3 -g -std=c++11 -ffast-math -funsafe-math-optimizations \
        -march=native -mtune=native -pedantic \
        -I/opt/local/include -I../src -I$(AURYNDIR)/src -I$(AURYNDIR)/src/auryn

LDFLAGS=$(BUILDDIR)/src/libauryn.a \
		-L/opt/local/lib -lboost_program_options-mt -lboost_serialization-mt -lboost_mpi-mt -lboost_system-mt -lboost_filesystem-mt

# debug
# BUILDDIR = $(AURYNDIR)/build/debug/
# CFLAGS= -ansi -Wall -pipe -g \
#         -march=native -mtune=native -pedantic \
#         -I$(AURYNDIR)/src -I$(AURYNDIR)/dev/src


OBJECTS=BiasIdentityConnection.o BiasConnection.o BurstPoissonGroup.o AdaptiveEBCPConnection.o EBCPConnection.o BCPConnection.o

SIMULATIONS=test_bias_ident_conn test_bias_conn test_burstpoisson test_burstpoisson_simple test_ficurve_burstpoisson test_adaptive_ebcp

.SECONDARY:

all: $(SIMULATIONS) 

test_%: test_%.o $(OBJECTS)
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
