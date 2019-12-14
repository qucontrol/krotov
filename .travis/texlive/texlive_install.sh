#!/usr/bin/env sh

# Originally from https://github.com/latex3/latex3

# This script is used for building LaTeX files using Travis
# A minimal current TL is installed adding only the packages that are
# required

# See if there is a cached version of TL available
if ! command -v texlua > /dev/null; then
  # Obtain TeX Live
  wget http://mirror.ctan.org/systems/texlive/tlnet/install-tl-unx.tar.gz
  tar -xzf install-tl-unx.tar.gz
  cd install-tl-20*

  # Install a minimal system
  ./install-tl --profile=../.travis/texlive/texlive.profile

  cd ..
fi

# Just including texlua so the cache check above works
tlmgr install luatex

# We specify the directory in which it is located texlive_packages
tlmgr install $(sed 's/\s*#.*//;/^\s*$/d' .travis/texlive/texlive_packages)

# Keep no backups (not required, simply makes cache bigger)
tlmgr option -- autobackup 0

# Update the TL install but add nothing new
tlmgr update --self --all --no-auto-install
