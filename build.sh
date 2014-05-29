#!/bin/sh

# If you're reading this, that means you've obtained this source code by cloning
# the repository instead of unpacking a distribution tarball.

# This script will change the modification date of certain files required by the
# GNU build system such that if your version of automake is less that 1.14 or
# your version of autoconf is less than 2.69, you won't have to upgrade just
# to build a fresh wnd-charm checkout.

touch aclocal.m4 Makefile.in configure config.h.in
./configure
make
