set(src_posthf
    src/pp_utilities.f90
    src/read_orbitals_from_file.f90
    src/twobody_hamiltonian.f90
    src/orbital_generators.f90
    src/onebody_hamiltonian.f90
    src/mp2_module.f90
    src/pyscf_driver.f90
    src/posthf_module.f90
   )
qe_enable_cuda_fortran("{src_posthf}")

qe_add_library(esh5_c src/esh5_interfaces_posthf.c)  
target_link_libraries(esh5_c
    PRIVATE
	qe_hdf5_c)

qe_add_library(qe_posthf ${src_posthf})
target_link_libraries(qe_posthf
    PRIVATE
        qe_openmp_fortran
        qe_pw
        qe_modules
        qe_fftx
        qe_upflib
        qe_lax
        qe_xclib
	esh5_c)

###########################################################


###########################################################
# posthf
###########################################################
set(src_posthf_x src/posthf.f90)
qe_add_executable(qe_posthf_exe ${src_posthf_x})
set_target_properties(qe_posthf_exe PROPERTIES OUTPUT_NAME posthf.x)
target_link_libraries(qe_posthf_exe
    PRIVATE
        qe_pw
        qe_modules
        qe_pp
        qe_fftx
        qe_upflib
        qe_xclib
	esh5_c
	qe_posthf)

###########################################################

qe_install_targets(
    # Libraries
    qe_posthf
    esh5_c
    # Executables
    qe_posthf_exe
    )
