##! /bin/bash
################################################################################
# Configure Cactus
################################################################################

# Pass configuration options to build script
echo "BEGIN MAKE_DEFINITION"
echo "EIDR_BUILD          = ${EIDR_BUILD}"
echo "EID_READ_EXTRA_LIB_DIRS = ${EIDR_EXTRA_LIB_DIRS}"
echo "EID_READ_EXTRA_LIBS     = ${EIDR_EXTRA_LIBS}"
echo "EID_READ_INSTALL_DIR    = ${EIDR_INSTALL_DIR}"
echo "END MAKE_DEFINITION"

# Set options
EID_READ_INC_DIRS="${EIDR_DIR}/include"
EID_READ_LIB_DIRS="${EIDR_DIR}/lib"
EID_READ_LIBS="${EIDR_EXTRA_LIBS}"

EID_READ_INC_DIRS="$(${CCTK_HOME}/lib/sbin/strip-incdirs.sh ${EIDR_INC_DIRS})"
EID_READ_LIB_DIRS="$(${CCTK_HOME}/lib/sbin/strip-libdirs.sh ${EIDR_LIB_DIRS})"

# Pass options to Cactus
echo "BEGIN MAKE_DEFINITION"
echo "EIDR_DIR      = ${EIDR_DIR}"
echo "EID_READ_INC_DIRS = ${EIDR_INC_DIRS}"
echo "EID_READ_LIB_DIRS = ${EIDR_LIB_DIRS}"
echo "EID_READ_LIBS     = ${EIDR_LIBS}"
echo "END MAKE_DEFINITION"

echo 'INCLUDE_DIRECTORY $(EID_READ_INC_DIRS)'
echo 'LIBRARY_DIRECTORY $(EID_READ_LIB_DIRS)'
echo 'LIBRARY           $(EID_READ_LIBS)'
