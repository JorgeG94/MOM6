# FetchFMS.cmake - FMS acquisition (find installed or fetch from git)
#
# MOM6_ACCESS=OFF (default): Fetch FMS from git (GFDL workflow)
# MOM6_ACCESS=ON:            Use find_package (ACCESS-NRI spack workflow)

include(FetchContent)

# Macro to find or fetch FMS
#
# Usage:
#   fetch_fms(URL <git_url> TAG <git_tag>)
#
macro(fetch_fms)
  cmake_parse_arguments(FMS "" "URL;TAG" "" ${ARGN})

  # Set defaults if not specified
  if(NOT FMS_URL)
    set(FMS_URL "https://github.com/NOAA-GFDL/FMS.git")
  endif()
  if(NOT FMS_TAG)
    set(FMS_TAG "ompoffload")
  endif()

  if(MOM6_ACCESS)
    # ACCESS-NRI workflow: use spack-installed FMS
    find_package(FMS REQUIRED CONFIG)
    message(STATUS "Found installed FMS: ${FMS_DIR}")

    # The installed FMS target doesn't export GPU link flags, but the FMS
    # static library contains GPU device code that requires the GPU runtime
    # at link time. Append these flags so they propagate to consumers.
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
      set_property(TARGET FMS::fms APPEND PROPERTY
        INTERFACE_LINK_OPTIONS -mp=gpu -acc=gpu -gpu=mem:separate)
    endif()

    # Create alias if needed (spack may export as FMS::fms)
    if(NOT TARGET FMS::fms_r8)
      if(TARGET FMS::fms)
        add_library(FMS::fms_r8 ALIAS FMS::fms)
      endif()
    endif()
  else()
    # GFDL workflow: fetch FMS from git
    message(STATUS "Fetching FMS from ${FMS_URL} (tag: ${FMS_TAG})")

    # Build FMS with its GPU offload flags, which also propagates the GPU
    # link options to MOM6 (only supported for NVHPC)
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
      set(FMS_OMPOFFLOAD ON CACHE BOOL "Enable GPU offload in FMS")
    endif()

    # Configure FMS build options
    set(OPENMP ${MOM6_OPENMP} CACHE BOOL "Enable OpenMP in FMS" FORCE)
    set(FMS_INSTALL OFF CACHE BOOL "" FORCE)

    # Fetch and build FMS
    FetchContent_Declare(fms
      GIT_REPOSITORY ${FMS_URL}
      GIT_TAG        ${FMS_TAG}
      GIT_SHALLOW    TRUE
      GIT_PROGRESS   TRUE
    )
    FetchContent_MakeAvailable(fms)

    # Create FMS::fms_r8 alias for compatibility
    if(NOT TARGET FMS::fms_r8)
      add_library(FMS::fms_r8 ALIAS fms)
    endif()

    message(STATUS "FMS fetched and configured successfully")
  endif()
endmacro()
