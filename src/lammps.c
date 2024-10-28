// LAMMPS 2001 - Molecular Dynamics Simulator
// Sandia National Laboratories, www.cs.sandia.gov/~sjplimp/lammps.html
// Steve Plimpton, sjplimp@sandia.gov
//
// Copyright (1998-2001) Sandia Corporation.  Under the terms of Contract
// DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
// certain rights in this software.  This software is distributed under 
// the GNU General Public License.
//
// See the README file in the top-level LAMMPS directory.

#include <stdio.h>

void initialize(void);
void input(int *iopflag);

int main (void)
{
	int iopflag = 0;
	initialize();
	input(&iopflag);
	return 0;
}
