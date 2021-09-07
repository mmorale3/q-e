# Makefile for postprocessing (PP)

include ../../make.inc

# location of needed modules and included files (if any)
MODFLAGS= $(BASEMOD_FLAGS) \
          $(MOD_FLAG)../../PW/src 

PPOBJS = \
posthf_module.o \
read_orbitals_from_file.o \
pp_utilities.o \
twobody_hamiltonian.o \
orbital_generators.o \
onebody_hamiltonian.o \
mp2_module.o \
pyscf_driver.o 

pyscf_driver_objs = \
posthf_module.o \
read_orbitals_from_file.o \
twobody_hamiltonian.o \
orbital_generators.o \
onebody_hamiltonian.o \
mp2_module.o \
pp_utilities.o

PWOBJS = ../../PW/src/libpw.a
QEMODS = ../../Modules/libqemod.a ../../upflib/libupf.a ../../KS_Solvers/libks_solvers.a \
         ../../FFTXlib/libqefft.a ../../LAXlib/libqela.a ../../UtilXlib/libutil.a \
         ../../dft-d3/libdftd3qe.a ../../GScratch/libgscratch.a 

MODULES = $(PWOBJS) $(QEOBJS)

TLDEPS= pwlibs

all : pw2posthf.x

libpp.a : $(PPOBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

pw2posthf.x : libpp.a pw2posthf.o $(MODULES) $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
                pw2posthf.o libpp.a ../../clib/clib.a $(MODULES) $(QEMODS) $(LIBOBJS) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../PP/src/$@ . )

pyscf_driver : libpp.a pyscf_driver.o $(pyscf_driver_objs) $(MODULES) $(QEMODS) $(LIBOBJS)
	f2py --f77exec=mpif90 --f90exec=mpif90 -c pyscf_driver.f90 -m pyscf_driver libpp.a ../../clib/clib.a \
                -I../../Modules/ $(MODFLAGS) $(MODULES) $(LIBOBJS) $(QELIBS2) \
		-L/usr/gapps/afqmc/libs/LASSEN/hdf5-1.10.5/lib /usr/gapps/afqmc/libs/LASSEN/hdf5-1.10.5/lib/libhdf5hl_fortran.a /usr/gapps/afqmc/libs/LASSEN/hdf5-1.10.5/lib/libhdf5_hl.a /usr/gapps/afqmc/libs/LASSEN/hdf5-1.10.5/lib/libhdf5_fortran.a /usr/gapps/afqmc/libs/LASSEN/hdf5-1.10.5/lib/libhdf5.a -lz -ldl -lm

clean :
	- /bin/rm -f *.x *.o *~ *_tmp.f90 *.d *.mod *.i *.L libpp.a

include make.depend
