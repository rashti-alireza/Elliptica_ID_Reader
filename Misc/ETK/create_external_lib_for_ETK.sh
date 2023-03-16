#!/usr/bin/env bash

##
## Creating an external library (thorn) for ETK (Cactus)
## Usage: (requires net access and git tool)
## $./me /path/to/cactus/dir
##
##
## NOTE: to update the Cactus branch of this repository:
## checkout the Cactus branch then issue:
## cp -vr /path/to/cactus/repos/ExternalLibraries-Elliptica_ID_Reader/* /path/to/Elliptica_ID_Reader/
##


if [[ $# -lt 1 ]]
then
	echo 'See the header of this script for the usage.'
	exit -1
fi

## conventions:
lib_version='1.0'
curdir=$(pwd)
cactus_root=$(realpath $1)
cactus_repo="${cactus_root}/repos"
cactus_arr="${cactus_root}/arrangements/ExternalLibraries"
lib_name='Elliptica_ID_Reader'
thorn="ExternalLibraries-Elliptica_ID_Reader"
thorn_path="${cactus_repo}/${thorn}"
git_repo="git@github.com:rashti-alireza/${lib_name}.git"
tarball_name="${lib_name}_${lib_version}"

## cleanup
cd ${cactus_repo}
rm -rf ${thorn}

## cp ExternalLibraries-Elliptica_ID_Reader from here to cactus
echo 'cp ...'
cd ${cactus_repo}
cp -r ${curdir}/${thorn} ${thorn}

## cloning
mkdir -p ${thorn_path}/dist
cd ${thorn_path}/dist
git clone --depth=1 ${git_repo}
cd ${thorn_path}/dist/${lib_name}
rm -rf .git
cd ${thorn_path}/dist
## note: it extract by this name
mv ${lib_name} ${tarball_name}
tar -zcf ${tarball_name}.tar.gz ${tarball_name}
rm -rf ${lib_name}

## set vars
cd ${thorn_path}/src

sed -i -E "s/^NAME=.+/NAME=${tarball_name}/g" build.sh

## set arrangments(these are compiled by cactus)
echo 'arranging ...'
cd ${cactus_arr}
rm -rf ${lib_name}
ln -s  ../../repos/ExternalLibraries-${lib_name} ${lib_name}


echo "Done! ($?)"





