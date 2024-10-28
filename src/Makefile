#!/usr/bin/make
#
# LMP				October 26, 2024
#
# source: Makefile
# author: @misael-diaz
#
# Synopsis:
# Defines the Makefile for building the program with GNU make.
#
# Copyright (c) 2024 Misael Diaz-Maldonado
# This file is released under the GNU General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#

include make-inc

all: $(GLOBAL_O)

$(GLOBAL_O):
	$(CC) $(CCOPT) -c $(GLOBAL_C) -o $(GLOBAL_O)

clean:
	rm -f *.o
