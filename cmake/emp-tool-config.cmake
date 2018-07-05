find_package(OpenSSL REQUIRED)
find_package(relic REQUIRED)
find_package(Boost REQUIRED COMPONENTS system)

# GMP
# https://raw.githubusercontent.com/stevedekorte/io/master/modules/FindGMP.cmake

# Try to find the GMP librairies
# GMP_FOUND - system has GMP lib
# GMP_INCLUDE_DIR - the GMP include directory
# GMP_LIBRARIES - Libraries needed to use GMP

if (GMP_INCLUDE_DIR AND GMP_LIBRARIES)
		# Already in cache, be silent
		set(GMP_FIND_QUIETLY TRUE)
endif (GMP_INCLUDE_DIR AND GMP_LIBRARIES)

find_path(GMP_INCLUDE_DIR NAMES gmp.h )
find_library(GMP_LIBRARIES NAMES gmp libgmp )
find_library(GMPXX_LIBRARIES NAMES gmpxx libgmpxx )
MESSAGE(STATUS "GMP libs: " ${GMP_LIBRARIES} " " ${GMPXX_LIBRARIES} )


mark_as_advanced(GMP_INCLUDE_DIR GMP_LIBRARIES)
# end of GMP

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
# End of MPFR




find_path(EMP-TOOL_INCLUDE_DIR NAMES emp-tool/emp-tool.h)
find_library(EMP-TOOL_LIBRARY NAMES emp-tool)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(EMP-TOOL DEFAULT_MSG EMP-TOOL_INCLUDE_DIR EMP-TOOL_LIBRARY)
add_definitions(-DEMP_CIRCUIT_PATH=${EMP-TOOL_INCLUDE_DIR}/emp-tool/circuits/files/)
if(EMP-TOOL_FOUND)
	set(EMP-TOOL_LIBRARIES ${EMP-TOOL_LIBRARY} ${RELIC_LIBRARIES} ${OPENSSL_LIBRARIES} ${Boost_LIBRARIES} ${GMP_LIBRARIES} ${MPFR_LIBRARIES})
	set(EMP-TOOL_INCLUDE_DIRS ${EMP-TOOL_INCLUDE_DIR} ${RELIC_INCLUDE_DIR} ${OPENSSL_INCLUDE_DIR} ${Boost_INCLUDE_DIRS} ${GMP_INCLUDE_DIR} ${MPFR_INCLUDE_DIR})
endif()
