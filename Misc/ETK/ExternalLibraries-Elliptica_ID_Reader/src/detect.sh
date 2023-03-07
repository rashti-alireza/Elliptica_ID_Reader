#! /bin/bash

################################################################################
# Prepare
################################################################################

# Set up shell
if [ "$(echo ${VERBOSE} | tr '[:upper:]' '[:lower:]')" = 'yes' ]; then
    set -x                      # Output commands
fi
set -e                          # Abort on errors



################################################################################
# Search
################################################################################

if [ -z "${Elliptica_ID_Reader_DIR}" ]; then
    echo "BEGIN MESSAGE"
    echo "Elliptica_ID_Reader selected, but Elliptica_ID_Reader_DIR not set."
    echo "END MESSAGE"
else
    echo "BEGIN MESSAGE"
    echo "Using Elliptica_ID_Reader in ${Elliptica_ID_Reader_DIR}"
    echo "END MESSAGE"
fi

THORN=Elliptica_ID_Reader


################################################################################
# Build
################################################################################

## here we set Elliptica_ID_Reader_DIR and INSTALL_DIR vars.
## NOTE: Elliptica_ID_Reader_INSTALL_DIR var managed by other defaults 
## setting of Cactus.
if [ -z "${Elliptica_ID_Reader_DIR}"                                                 \
     -o "$(echo "${Elliptica_ID_Reader_DIR}" | tr '[a-z]' '[A-Z]')" = 'BUILD' ]
then
    echo "BEGIN MESSAGE"
    echo "Using bundled Elliptica_ID_Reader..."
    echo "END MESSAGE"
    
    # Check for required tools. Do this here so that we don't require
    # them when using the system library.
    if [ "x$TAR" = x ] ; then
        echo 'BEGIN ERROR'
        echo 'Could not find tar command.'
        echo 'Please make sure that the (GNU) tar command is present,'
        echo 'and that the TAR variable is set to its location.'
        echo 'END ERROR'
        exit 1
    fi
    if [ "x$PATCH" = x ] ; then
        echo 'BEGIN ERROR'
        echo 'Could not find patch command.'
        echo 'Please make sure that the patch command is present,'
        echo 'and that the PATCH variable is set to its location.'
        echo 'END ERROR'
        exit 1
    fi

    # Set locations
    SRCDIR="$(dirname $0)"
    BUILD_DIR=${SCRATCH_BUILD}/build/${THORN}
    if [ -z "${Elliptica_ID_Reader_INSTALL_DIR}" ]; then
        INSTALL_DIR=${SCRATCH_BUILD}/external/${THORN}
    else
        echo "BEGIN MESSAGE"
        echo "Installing Elliptica_ID_Reader into ${Elliptica_ID_Reader_INSTALL_DIR}"
        echo "END MESSAGE"
        INSTALL_DIR=${Elliptica_ID_Reader_INSTALL_DIR}
    fi
    Elliptica_ID_Reader_BUILD=1
    Elliptica_ID_Reader_DIR=${INSTALL_DIR}
else
    Elliptica_ID_Reader_BUILD=
    DONE_FILE=${SCRATCH_BUILD}/done/${THORN}
    if [ ! -e ${DONE_FILE} ]; then
        mkdir ${SCRATCH_BUILD}/done 2> /dev/null || true
        date > ${DONE_FILE}
    fi
fi



################################################################################
# Configure Cactus
################################################################################

# Pass configuration options to build script
echo "BEGIN MAKE_DEFINITION"
echo "Elliptica_ID_Reader_BUILD          = ${Elliptica_ID_Reader_BUILD}"
echo "Elliptica_ID_Reader_INSTALL_DIR    = ${Elliptica_ID_Reader_INSTALL_DIR}"
echo "END MAKE_DEFINITION"

# Set options
Elliptica_ID_Reader_INC_DIRS="${Elliptica_ID_Reader_DIR}/include"
Elliptica_ID_Reader_LIB_DIRS="${Elliptica_ID_Reader_DIR}/lib"
Elliptica_ID_Reader_LIBS="elliptica_id_reader"
## rename things for cactus
Elliptica_ID_Reader_INC_DIRS="$(${CCTK_HOME}/lib/sbin/strip-incdirs.sh ${Elliptica_ID_Reader_INC_DIRS})"
Elliptica_ID_Reader_LIB_DIRS="$(${CCTK_HOME}/lib/sbin/strip-libdirs.sh ${Elliptica_ID_Reader_LIB_DIRS})"

# Pass options to Cactus
echo "BEGIN MAKE_DEFINITION"
echo "Elliptica_ID_Reader_DIR      = ${Elliptica_ID_Reader_DIR}"
echo "Elliptica_ID_Reader_INC_DIRS = ${Elliptica_ID_Reader_INC_DIRS}"
echo "Elliptica_ID_Reader_LIB_DIRS = ${Elliptica_ID_Reader_LIB_DIRS}"
echo "Elliptica_ID_Reader_LIBS     = ${Elliptica_ID_Reader_LIBS}"
echo "END MAKE_DEFINITION"

echo 'INCLUDE_DIRECTORY $(Elliptica_ID_Reader_INC_DIRS)'
echo 'LIBRARY_DIRECTORY $(Elliptica_ID_Reader_LIB_DIRS)'
echo 'LIBRARY           $(Elliptica_ID_Reader_LIBS)'

