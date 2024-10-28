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

all: $(LMP)

$(GLOBAL_O): $(GLOBAL_C)
	$(CC) $(CCOPT) -c $(GLOBAL_C) -o $(GLOBAL_O)
$(INITIALIZE_O): $(INITIALIZE_C)
	$(CC) $(CCOPT) -c $(INITIALIZE_C) -o $(INITIALIZE_O)
$(INPUT_O): $(INPUT_C)
	$(CC) $(CCOPT) -c $(INPUT_C) -o $(INPUT_O)
$(LAMMPS_O): $(LAMMPS_C)
	$(CC) $(CCOPT) -c $(LAMMPS_C) -o $(LAMMPS_O)
$(LMP): $(GLOBAL_O) $(INITIALIZE_O) $(INPUT_O) $(LAMMPS_O)
	$(CC) $(CCOPT) $(GLOBAL_O) $(INITIALIZE_O) $(INPUT_O) $(LAMMPS_O) -o $(LMP)

clean:
	rm -f *.o lmp
