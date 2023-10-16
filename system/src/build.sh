##! /bin/bash
################################################################################
# Configure Cactus
################################################################################

# Pass configuration options to build script
echo "BEGIN MAKE_DEFINITION"
echo "EIDR_BUILD          = ${EIDR_BUILD}"
echo "EIDRBHNS_EXTRA_LIB_DIRS = ${EIDR_EXTRA_LIB_DIRS}"
echo "EIDRBHNS_EXTRA_LIBS     = ${EIDR_EXTRA_LIBS}"
echo "EIDRBHNS_INSTALL_DIR    = ${EIDR_INSTALL_DIR}"
echo "END MAKE_DEFINITION"

# Set options
EIDRBHNS_INC_DIRS="${EIDR_DIR}/include"
EIDRBHNS_LIB_DIRS="${EIDR_DIR}/lib"
EIDRBHNS_LIBS="${EIDR_EXTRA_LIBS}"

EIDRBHNS_INC_DIRS="$(${CCTK_HOME}/lib/sbin/strip-incdirs.sh ${EIDR_INC_DIRS})"
EIDRBHNS_LIB_DIRS="$(${CCTK_HOME}/lib/sbin/strip-libdirs.sh ${EIDR_LIB_DIRS})"

# Pass options to Cactus
echo "BEGIN MAKE_DEFINITION"
echo "EIDR_DIR      = ${EIDR_DIR}"
echo "EIDRBHNS_INC_DIRS = ${EIDR_INC_DIRS}"
echo "EIDRBHNS_LIB_DIRS = ${EIDR_LIB_DIRS}"
echo "EIDRBHNS_LIBS     = ${EIDR_LIBS}"
echo "END MAKE_DEFINITION"

echo 'INCLUDE_DIRECTORY $(EIDRBHNS_INC_DIRS)'
echo 'LIBRARY_DIRECTORY $(EIDRBHNS_LIB_DIRS)'
echo 'LIBRARY           $(EIDRBHNS_LIBS)'
