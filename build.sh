#!/bin/sh

touch aclocal.m4 Makefile.in configure config.h.in
./configure
make
