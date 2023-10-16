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
# Build
################################################################################

# Set locations
THORN='Elliptica_ID_Reader'
NAME='Elliptica_ID_Reader'

SRCDIR="$(dirname $0)"
BUILD_DIR=${SCRATCH_BUILD}/build/${THORN}
DONE_FILE=${SCRATCH_BUILD}/done/${THORN}
if [ -z "${Elliptica_ID_Reader_INSTALL_DIR}" ]; then
    INSTALL_DIR=${SCRATCH_BUILD}/external/${THORN}
    echo "BEGIN MESSAGE"
    echo "Installing Elliptica_ID_Reader in the default ${INSTALL_DIR} directory."
    echo "END MESSAGE"
else
    echo "BEGIN MESSAGE"
    echo "Installing Elliptica_ID_Reader into ${Elliptica_ID_Reader_INSTALL_DIR}"
    echo "END MESSAGE"
    INSTALL_DIR==${Elliptica_ID_Reader_INSTALL_DIR}
fi
Elliptica_ID_Reader_DIR=${INSTALL_DIR}

# Set up environment
unset LIBS
if echo '' ${ARFLAGS} | grep 64 > /dev/null 2>&1; then
    export OBJECT_MODE=64
fi

echo "Elliptica_ID_Reader: Preparing directory structure..."
cd ${SCRATCH_BUILD}
mkdir build external done 2> /dev/null || true
rm -rf ${BUILD_DIR} ${INSTALL_DIR}
mkdir ${BUILD_DIR} ${INSTALL_DIR}

#echo "Elliptica_ID_Reader: Unpacking archive..."
pushd ${BUILD_DIR}
#${TAR?} xzf ${SRCDIR}/../dist/${NAME}.tar.gz
cp -r "${SRCDIR}/../dist/${NAME}" .
cd ${NAME}

#echo "Elliptica_ID_Reader: Configuring..."
## Elliptica_ID_Reader doesn't need any configuration script

echo "Elliptica_ID_Reader: Building..."
## NOTE: these vars coming from Cactus
(
##    exec >/dev/null # redirect stdout to null to avoid makefile prints
    ${MAKE} CC="${CC}" CFLAGS="${CFLAGS}" AR="${AR}" LIB_TYPE="static"
)
if (( $? )); then
    echo 'BEGIN ERROR'
    echo 'Error while building Elliptica_ID_Reader. Aborting.'
    echo 'END ERROR'
    exit 1
fi

mv -v ./include ${INSTALL_DIR}/
mv -v ./lib ${INSTALL_DIR}/

popd

echo "Elliptica_ID_Reader: Cleaning up..."
rm -rf ${BUILD_DIR}

date > ${DONE_FILE}
echo "Elliptica_ID_Reader: Done."
