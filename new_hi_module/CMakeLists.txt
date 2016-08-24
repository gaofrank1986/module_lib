# CMake project file for NEMOH 
 
cmake_minimum_required (VERSION 2.8) 
project (dutwav Fortran) 
enable_language (Fortran) 
 
if(NOT CMAKE_BUILD_TYPE) 
	SET(CMAKE_BUILD_TYPE Release CACHE STRING "default to debug" FORCE) 
endif(NOT CMAKE_BUILD_TYPE) 
 
SET(dutwav_LIB_TYPE "SHARED") 
 
IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin") 
	set(CMAKE_MACOSX_RPATH ON) 
	SET(CMAKE_SKIP_RPATH  TRUE) 
	SET(CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES}) 
	SET(BUILD_SHARED_LIBRARIES ON) 
	SET(CMAKE_EXE_LINKER_FLAGS "-static -static-libgfortran -static-libgcc") 
	SET(CMAKE_SHARED_LINKER_FLAGS "-static -static-libgfortran -static-libgcc ") 
	SET(CMAKE_MODULE_LINKER_FLAGS "-static -static-libgfortran -static-libgcc") 
    #SET(NEMOH_LIB_TYPE "STATIC") 
 
ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin") 
 
if(WIN32) 
 
	if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel") 
		set(CMAKE_Fortran_FLAGS "/stand /Qopenmp /Qmkl ") 
		set(CMAKE_Fortran_FLAGS_DEBUG "/stand /Qopenmp /Qmkl /fp:source /Od /traceback /check:all /warn:all /debug:all ") 
		set(CMAKE_Fortran_FLAGS_RELEASE "/stand /Qopenmp /Qmkl ") 
	 
		# You can uncomment the following lines for further optimisation but doing so results will not be completely exact with original version. 
		#set(CMAKE_Fortran_FLAGS "/stand /Qopenmp /Qmkl /fast") 
		 
	elseif ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU") 
		set(CMAKE_Fortran_FLAGS "-W -Wall -std=f2008 -fopenmp -fmax-errors=0  -O3 -w -g ") 
		set(CMAKE_Fortran_FLAGS_DEBUG "-W -Wall -std=f2008 -fopenmp   -fmax-errors=0 -O3 -w -g ") 
		set(CMAKE_Fortran_FLAGS_RELEASE "-W -Wall -std=f2008 -fopenmp   -fmax-errors=0 -w") 
	else() 
	 
		message(SEND_ERROR "Only gfortran and Intel ifort compiler are currently supported.") 
	 
	endif() 
endif (WIN32) 
 
if (UNIX) 
	if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel") 
		set(CMAKE_Fortran_FLAGS "-stand   -heap-arrays 1024") 
		#set(CMAKE_Fortran_FLAGS "-stand -openmp -mkl -heap-arrays 1024") 
		set(CMAKE_Fortran_FLAGS_DEBUG "-stand   -traceback -g -fp-model source -check all -warn all -debug all -heap-arrays 1024") 
		set(CMAKE_Fortran_FLAGS_RELEASE "-stand  -heap-arrays 1024") 
	 
		# You can uncomment the following lines for further optimisation but doing so results will not be completely exact with original version. 
		#set(CMAKE_Fortran_FLAGS "-stand -openmp -mkl -fast") 
		 
	elseif ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU") 
		set(CMAKE_Fortran_FLAGS "-W -Wall -std=f2008 -fopenmp   -fmax-errors=0  -O3 -w -g ") 
		set(CMAKE_Fortran_FLAGS_DEBUG "-W -Wall -std=f2008 -fopenmp   -fmax-errors=0  -O3 -w -g   ") 
		set(CMAKE_Fortran_FLAGS_RELEASE "-W -Wall -std=f2008 -fopenmp   -fmax-errors=0 -O3 -w  ") 
	 
	else() 
	 
		message(SEND_ERROR "Only gfortran and Intel ifort compiler are currently supported.") 
	endif() 
 
	find_package(BLAS REQUIRED) 
	find_package(LAPACK REQUIRED) 
endif (UNIX) 
 
 
 
 
 
set(CMAKE_EXE_LINKER_FLAGS "") 

 
if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU") 
 
	SET_TARGET_PROPERTIES(nemoh PROPERTIES COMPILE_FLAGS "-fopenmp" LINK_FLAGS "-fopenmp") 
	SET_TARGET_PROPERTIES(nemoh PROPERTIES COMPILE_FLAGS "-fopenmp -lblas -llapack" LINK_FLAGS "-fopenmp -lblas -llapack") 
	target_link_libraries(nemoh ${BLAS_LIBRARIES}) 
	target_link_libraries(nemoh ${LAPACK_LIBRARIES}) 
	target_link_libraries(nemoh ${OpenMP_LIBRARIES}) 
endif()


ADD_SUBDIRECTORY(src)

add_executable(test "./test.f90" )
target_link_libraries(test hi2_mod)
