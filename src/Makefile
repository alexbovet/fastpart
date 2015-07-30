PROG =	fastpart

SRCS =  main.f90 basic_mod.f90 nrutil.f90 newrun.f90 restart.f90 \
        auxval.f90 initial.f90 resume.f90 start.f90 diagnose.f90 boris.f90\
        stepon.f90 tesend.f90 endrun.f90 chkrst.f90 mv2bk.f90 memory.f90 interpolate.f90 \
        injection_mod.f90 diagnostic_mod.f90 diagnose_1d.f90 numerics_mod.f90 tracers_mod.f90 \
        fields_mod.f90 cart_to_cyl.f90 cyl_to_cart.f90 cyl_to_tilt.f90 tilt_to_cyl.f90 \
	Bfieldfind.f90 Efieldfind.f90 sinegrid.f90 GBSload.f90 centdiff.f90 rk4_3d.f90 random_mod.f90 \
        constants.f90

SRCS_C =

OBJS = ${SRCS:.f90=.o} ${SRCS_F90:.F90=.o} ${SRCS_C:.c=.o}

LIBS = -lfutils \
      -lhdf5_fortran -lhdf5 -lz \
      -liomp5 -lpthread \
      -lpthread

CC = icc
CFLAGS =
FC = mpif90

# OPTIMIZED 
#FFLAGS = -O3 -xSSE4.1 -I/usr/local/crpp/include/O
# OPTIMIZED for pleiades2
FFLAGS = -O3  -I${HOME}/include  

# OPTIMIZED AND PROFILE OUTPUT
#FFLAGS = -p -O3 -xSSE4.1 -I/usr/local/crpp/include/O
# FOR DEBUGGING AND TRACEBACK
#FFLAGS = -g -traceback -I/usr/local/crpp/include/g

F90 = $(FC)
F90FLAGS = $(FFLAGS)

# NO PROFILE OUTPUT
#LDFLAGS = -Wl,-rpath=$(MKL) -L/usr/local/crpp/lib/O -L${HDF5}/lib -L$(MKL)
# NO PROFILE OUTPUT for pleiades2
#LDFLAGS = -i-dynamic -Wl,-rpath=$(MKL) -L/usr/local/crpp/lib/O -L${HDF5}/lib -L$(MKL)
# PROFILE OUTPUT
LDFLAGS = -p -i-dynamic -Wl,-rpath=$(MKL) -L${HOME}/lib/O -L${HDF5}/lib 

all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

tags:
	etags *.f90

clean:
	rm -f $(OBJS) *.mod

distclean: clean
	rm -f $(PROG) *~ a.out *.o

.SUFFIXES: $(SUFFIXES) .f90

.f90.o:
	$(F90) $(F90FLAGS) -c $<

# some object files are linked in libraries above
main.o:		   basic_mod.o fields_mod.o
auxval.o:          basic_mod.o numerics_mod.o fields_mod.o tracers_mod.o
basic_mod.o:	   constants.o
boris.o:           basic_mod.o tracers_mod.o fields_mod.o numerics_mod.o Bfieldfind.o
cart_to_cyl.o:	   constants.o fields_mod.o diagnostic_mod.o tracers_mod.o injection_mod.o 
cyl_to_cart.o:     fields_mod.o diagnostic_mod.o tracers_mod.o injection_mod.o 
cyl_to_tilt.o:     fields_mod.o diagnostic_mod.o tracers_mod.o injection_mod.o 
tilt_to_cyl.o:     fields_mod.o diagnostic_mod.o tracers_mod.o injection_mod.o 
chkrst.o:          basic_mod.o
diagnose.o:        basic_mod.o diagnostic_mod.o tracers_mod.o fields_mod.o injection_mod.o numerics_mod.o
diagnose_1d.o:     basic_mod.o diagnostic_mod.o tracers_mod.o injection_mod.o fields_mod.o
endrun.o:          basic_mod.o
Bfieldfind.o:      basic_mod.o diagnostic_mod.o tracers_mod.o fields_mod.o injection_mod.o constants.o interpolate.o
Efieldfind.o:      basic_mod.o diagnostic_mod.o tracers_mod.o fields_mod.o injection_mod.o random_mod.o constants.o interpolate.o
sinegrid.o:	   basic_mod.o fields_mod.o constants.o
GBSload.o:	   basic_mod.o fields_mod.o constants.o
centdiff.o:        basic_mod.o fields_mod.o constants.o interpolate.o
rk4_3d.o:	   basic_mod.o tracers_mod.o fields_mod.o numerics_mod.o
initial.o:         constants.o basic_mod.o injection_mod.o tracers_mod.o random_mod.o
memory.o:          injection_mod.o tracers_mod.o numerics_mod.o fields_mod.o
stepon.o:          basic_mod.o numerics_mod.o tracers_mod.o fields_mod.o interpolate.o random_mod.o
tesend.o:          basic_mod.o
newrun.o:          basic_mod.o fields_mod.o diagnostic_mod.o injection_mod.o numerics_mod.o tracers_mod.o constants.o

