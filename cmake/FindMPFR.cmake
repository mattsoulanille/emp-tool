# https://raw.githubusercontent.com/stevedekorte/io/master/modules/FindGMP.cmake

# Try to find the MPFR librairies
# MPFR_FOUND - system has MPFR lib
# MPFR_INCLUDE_DIR - the MPFR include directory
# MPFR_LIBRARIES - Libraries needed to use MPFR

if (MPFR_INCLUDE_DIR AND MPFR_LIBRARIES)
		# Already in cache, be silent
		set(MPFR_FIND_QUIETLY TRUE)
endif (MPFR_INCLUDE_DIR AND MPFR_LIBRARIES)

find_path(MPFR_INCLUDE_DIR NAMES mpfr.h )
find_library(MPFR_LIBRARIES NAMES mpfr libmpfr )
find_library(MPFRXX_LIBRARIES NAMES mpfrxx libmpfrxx )
MESSAGE(STATUS "MPFR libs: " ${MPFR_LIBRARIES} " " ${MPFRXX_LIBRARIES} )

include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(MPFR DEFAULT_MSG MPFR_INCLUDE_DIR MPFR_LIBRARIES)

mark_as_advanced(MPFR_INCLUDE_DIR MPFR_LIBRARIES)
