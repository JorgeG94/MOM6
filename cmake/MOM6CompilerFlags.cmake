# MOM6 Compiler Flags
#
# This module defines compiler flags for MOM6 based on the compiler vendor.
# Flags are aligned with FMS cmake flags where possible, with MOM6-specific
# additions (big-endian conversion, backtrace).
#
# Categories:
#   FLAGS_STANDARD  - Always applied (real8, line length, etc.)
#   FLAGS_DEBUG     - Debug builds
#   FLAGS_RELEASE   - Optimized builds

# Helper function to join list elements into a string
function(list_to_string LIST_VAR OUTPUT_VAR)
  string(REPLACE ";" " " _result "${${LIST_VAR}}")
  set(${OUTPUT_VAR} "${_result}" PARENT_SCOPE)
endfunction()

# GNU Fortran (gfortran)
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(FLAGS_STANDARD
    -ffree-line-length-none
    -fconvert=big-endian
    -fdefault-real-8
    -fdefault-double-8
    -fcray-pointer
    -fbacktrace
  )

  set(FLAGS_DEBUG
    -O0
    -g
    -fcheck=bounds
  )

  set(FLAGS_RELEASE
    -O2
  )

  # GCC 10+ requires -fallow-argument-mismatch for legacy MPI interfaces
  if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL "10.0")
    list(APPEND FLAGS_STANDARD -fallow-argument-mismatch -fallow-invalid-boz)
  endif()

# Intel Classic Fortran (ifort)
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  set(FLAGS_STANDARD
    -fpp
    -fno-alias
    -safe-cray-ptr
    -ftz
    -assume byterecl
    -convert big_endian
    -traceback
    -r8
  )

  set(FLAGS_DEBUG
    -O0
    -g
    -check bounds
  )

  set(FLAGS_RELEASE
    -O2
    -debug minimal
  )

# Intel LLVM Fortran (ifx)
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
  set(FLAGS_STANDARD
    -fpp
    -fno-alias
    -safe-cray-ptr
    -ftz
    -assume byterecl
    -convert big_endian
    -traceback
    -r8
  )

  set(FLAGS_DEBUG
    -O0
    -g
    -check bounds
  )

  set(FLAGS_RELEASE
    -O2
    -debug minimal
  )

# NVIDIA HPC Fortran (nvfortran)
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
  set(FLAGS_STANDARD
    -O2
    -Mnovect
    -mp=gpu
    -stdpar=gpu
    -Mnofma
    -Minfo=accel
    -acc=gpu
    -gpu=mem:separate
    -Minline=exclude:${CMAKE_SOURCE_DIR}/src/framework
    # Currently required by NVIDIA GPU in MOM_continuity_PPM.F90
    -Minline=name:flux_elem,name:flux_elem_OBC,name:ratio_max
    -r8
  )

  # Optional GPU compute capability, e.g. -DMOM6_GPU_ARCH=cc90
  if(MOM6_GPU_ARCH)
    if(NOT MOM6_GPU_ARCH MATCHES "^cc[a-z0-9]+(,cc[a-z0-9]+)*$")
      message(FATAL_ERROR "MOM6_GPU_ARCH='${MOM6_GPU_ARCH}' is invalid; expected e.g. cc90, cc80,cc90 or ccnative")
    endif()
    list(APPEND FLAGS_STANDARD -gpu=${MOM6_GPU_ARCH})
    add_link_options(-gpu=${MOM6_GPU_ARCH})
    message(STATUS "MOM6 GPU arch: ${MOM6_GPU_ARCH}")
  endif()

  add_compile_definitions(
    HAVE_MPI=1
    HAVE_FC_DO_CONCURRENT_LOCALITY=1
  )

  set(FLAGS_DEBUG
    -O0
    -g
    -Mbounds
  )

  set(FLAGS_RELEASE
    -O2
  )

# Cray Fortran (crayftn)
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Cray")
  set(FLAGS_STANDARD
    -f free
    -h byteswapio
    -s real64
    -eF
  )

  set(FLAGS_DEBUG
    -O0
    -g
    -R bps
  )

  set(FLAGS_RELEASE
    -O2
  )

else()
  message(WARNING "Unknown Fortran compiler: ${CMAKE_Fortran_COMPILER_ID}. Using default flags.")
  set(FLAGS_STANDARD "")
  set(FLAGS_DEBUG "-O0 -g")
  set(FLAGS_RELEASE "-O2")
endif()

# Convert lists to space-separated strings
list_to_string(FLAGS_STANDARD FLAGS_STANDARD_STR)
list_to_string(FLAGS_DEBUG FLAGS_DEBUG_STR)
list_to_string(FLAGS_RELEASE FLAGS_RELEASE_STR)

# Apply flags globally to CMAKE_Fortran_FLAGS.  This is a normal (non-cache)
# variable so that the standard flags are not re-appended to the cached user
# value on every reconfigure.
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FLAGS_STANDARD_STR}")
set(CMAKE_Fortran_FLAGS_DEBUG "${FLAGS_DEBUG_STR}"
    CACHE STRING "Fortran compiler flags for Debug builds" FORCE)
set(CMAKE_Fortran_FLAGS_RELEASE "${FLAGS_RELEASE_STR}"
    CACHE STRING "Fortran compiler flags for Release builds" FORCE)
set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "${FLAGS_RELEASE_STR} -g"
    CACHE STRING "Fortran compiler flags for RelWithDebInfo builds" FORCE)

# Report the flags being used
message(STATUS "MOM6 Fortran Compiler: ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION}")
message(STATUS "MOM6 Fortran Flags: ${FLAGS_STANDARD_STR}")
