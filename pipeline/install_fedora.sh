#!/bin/bash
# An small script for libraries intallation
#
# Google's Shell Style Guide

# sudo dnf install redhat-rpm-config


function upgrade_system {
  sudo dnf install python-virtualenv
  sudo dnf install python-devel
}


function install_virtualenv {
  # Install virtualenv for deploy a new enviroment
  virtualenv --python=/usr/bin/python2.7  $installation_dir/.venv
  source $installation_dir/.venv/bin/activate
}


# Update python modules
function update_pip {
  pip install --upgrade pip
  pip install -r requirements.txt
}


function install_atlas {
  atlas_url="https://downloads.sourceforge.net/project/math-atlas/Stable/3.10.3/atlas3.10.3.tar.bz2?r=&ts=1506698217&use_mirror=10gbps-io"
  lapack_url="http://www.netlib.org/lapack/lapack-3.7.1.tgz"

  wget -O atlas.tar.bz2 $atlas_url
  wget -O lapack.tgz $lapack_url

  tar -xf atlas.tar.bz2
  rm atlas.tar.bz2

  cd ATLAS
  if [ ! -d DONE/ ]; then
    mkdir DONE
  fi
  cd DONE/

  mkdir /home/user/Work/Projects/pipeline/.local/ATLAS
  ../configure --shared -Fa alg -fPIC\
  --with-netlib-lapack-tarfile=../../lapack.tgz\
  --prefix=/home/user/Work/Projects/pipeline/.local/ATLAS

  mkdir /home/user/Work/Projects/pipeline/.local/ATLAS/lib

  sed -i -e 's/"-rpath-link $(LIBINSTdir)"/-rpath-link $(LIBINSTdir)/g' Makefile

  # Get into lib folder and make shared libraries
  cd lib
  # Change reference to fortran lib
  sed -i -e 's/4.8.5/4.8.2/g' Make.inc
  # Change ""
  sed -i -e 's/"-rpath-link $(LIBINSTdir)"/-rpath-link $(LIBINSTdir)/g' Makefile

  # make shared
  cd ../

  # Compile
  make
  make install
}


function install_cdsclient {
  cdsclient_url="http://cdsarc.u-strasbg.fr/ftp/pub/sw/cdsclient.tar.gz"
  wget -O cdsclient.tar.gz $cdsclient_url

  if [ ! -d cdsclient/ ]; then
    mkdir cdsclient/
  fi

  tar -xzf cdsclient.tar.gz -C cdsclient --strip-components=1
  rm cdsclient.tar.gz  # remove old file
  cd cdsclient

  # Create a dir for cdsclient installation
  cdsclient_dir="/home/user/Work/Projects/pipeline/.local/cdsclient"

  if [ ! -d $cdsclient_dir ]; then
    mkdir $cdsclient_dir
  fi

  # Configure
  ./configure -prefix=/home/user/Work/Projects/pipeline/.local/cdsclient

  # Compile them
  make
  # Perfom an installation to local folder
  make install
  # Remove old directory
  rm -r cdsclient
}


function install_sextractor {
  sextractor_url="https://www.astromatic.net/download/sextractor/sextractor-2.19.5.tar.gz"
  wget -O sextractor.tar.gz $sextractor_url

  if [ ! -d sextractor/ ]; then
    mkdir sextractor/
  fi

  tar -xzf sextractor.tar.gz -C sextractor/ --strip-components=1
  rm sextractor.tar.gz  # Remove old file
  cd sextractor

  # Configure
  ./configure --with-atlas-incdir=$1 --with-atlas-libdir=$2 --prefix=$3

  # Compile Sextractor
  make
  # Install Sextractor in local directories
  make install
  # Remove old directory
  rm -r sextractor
}


function install_scamp {
  scamp_url="https://www.astromatic.net/download/scamp/scamp-2.0.4.tar.gz"
  wget -O scamp.tar.gz $scamp_url

  if [ ! -d scamp/ ]; then
    mkdir scamp/
  fi

  tar -xzf scamp.tar.gz -C scamp/ --strip-components=1
  rm scamp.tar.gz  # Remove old file
  cd scamp

  # Configure
  cdsclient_bin_dir="/home/user/Work/Projects/pipeline/.local/cdsclient/bin"
  ./configure --with-atlas-incdir=$1 --with-atlas-libdir=$2\
  --with-cdsclient-dir=$cdsclient_bin_dir --prefix=$3

  # Compile Scamp
  make
  # Install Scamp in local directories
  make install
  # Remove old directory
  rm -r scamp
}


function update_enviroment {
  cp /media/sf_Euclid-tests/pipeline/.bash_profile ~/.bash_profile
  source ~/.bash_profile
}


# TODO improve format!
function copy_files {
  cp -r /media/sf_Euclid-tests/pipeline/* /home/user/Work/Projects/pipeline/
  cp -r /media/sf_Euclid-tests/pipeline/.settings.ini ~/Work/Projects/pipeline/.settings.ini
}


function main {
  tmp_dir="/home/user/Work/Projects/pipeline/tmp/"
  local_dir="/home/user/Work/Projects/pipeline/.local/"
  installation_dir="/media/sf_Euclid-tests/pipeline/"
  atlas_include_dir="/home/user/Work/Projects/pipeline/.local/ATLAS/include"
  atlas_lib_dir="/home/user/Work/Projects/pipeline/.local/ATLAS/lib"

  cd $installation_dir

  upgrade_system
  install_virtualenv
  update_pip

  # Checking directories
  if [ ! -d "$tmp_dir" ]; then
    mkdir $tmp_dir
  fi

  if [ ! -d "$local_dir" ]; then
    mkdir $local_dir
  fi

  # Install scamp from scratch
  # Compile ATLAS/Lapack library
  cd $tmp_dir

  install_atlas

  cd ../../
  rm -rf a*

  install_cdsclient
  cd ../

  install_sextractor $atlas_include_dir $atlas_lib_dir $local_dir
  cd ../

  install_scamp $atlas_include_dir $atlas_lib_dir $local_dir
  cd ../

  update_enviroment
  copy_files
}


if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
  main "$@"
fi
