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

// -------------------------------------------------------------------------
// read all info from text input file one line at a time
// -------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

static char *str;
static size_t sz_str;

void input (int *iopflag)
{
	FILE *f = fopen("data.txt", "r");
	if (!f) {
		fprintf(stderr, "input: %s\n", strerror(errno));
		exit(EXIT_FAILURE);
	}
	*iopflag = 0;
	int nlines = 0;
	while (1) {
		errno = 0;
		ssize_t rc = getline(&str, &sz_str, f);
		if (-1 == rc) {
			if (errno) {
				if (str) {
					free(str);
				}
				fclose(f);
				fprintf(stderr, "input: %s\n", strerror(errno));
				exit(EXIT_FAILURE);
			}
			*iopflag = EOF;
			break;
		}
		fprintf(stdout, "%s", str);
		++nlines;
	}
	fprintf(stdout, "nlines: %d\n", nlines);
	free(str);
	sz_str = 0;
	fclose(f);
}
