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

#define VERBOSE 1

extern int firstflag;

extern int units;
extern double boltz;
extern double dtfactor;
extern double pfactor;
extern double efactor;
extern double skin;
extern double cutlj;
extern double cutcoul;

extern int idimension;
extern double dt;

extern int nonstyle;
extern int offsetflag;
extern int mixflag;
extern int mixstyle;

static char *str;
static size_t sz_str;

static void in_units (char *string)
{
	char const *str = string;
	if (strstr(str, "real")) {
		units = 0;
		boltz = 0.001987191;
		dtfactor = 48.88821;
		pfactor = 68589.796;
		efactor = 332.0636;
		skin = 2.0;
		cutlj = 10.0;
		cutcoul = 10.0;
		fprintf(stdout, "Units real\n");
	} else if (strstr(str, "lj")) {
		units = 1;
		boltz = 1;
		dtfactor = 1;
		pfactor = 1;
		efactor = 1;
		skin = 0.3;
		cutlj = 2.5;
		cutcoul = 2.5;
		fprintf(stdout, "input: Units lj\n");
	} else {
		free(string);
		fprintf(stderr, "input: bad units param\n");
		exit(EXIT_FAILURE);
	}
}

static void in_dimension (char *string)
{
	char const *str = string;
	char const *s = strstr(str, "dimension");
	idimension = atoi(&s[9]);
	if (2 != idimension && 3 != idimension) {
		free(string);
		fprintf(stderr, "input: bad dimension param\n");
		exit(EXIT_FAILURE);
	}
	fprintf(stdout, "input: Dimension %d\n", idimension);
}

static void in_timestep (char *string)
{
	char const *str = string;
	char const *prm = "timestep";
	char const *s = strstr(str, prm);
	dt = atof(&s[strlen(prm)]);
	if (0 >= dt) {
		free(string);
		fprintf(stderr, "input: bad timestep param\n");
		exit(EXIT_FAILURE);
	}
	dt /= dtfactor;
	fprintf(stdout, "input: Timestep %le\n", dt);
}

static void in_nonbond_style (char *string)
{
	char const *str = string;
	char const *prm = "nonbond style";
	char const *none = "none";
	char const *ljcutoff = "lj/cutoff";
	char const *s = strstr(str, prm);
	if (strstr(s, none)) {
		nonstyle = 0;
	} else if (strstr(s, ljcutoff)) {

		s = strstr(s, ljcutoff);
		s = &s[strlen(ljcutoff)];

		nonstyle = 1;

		char *endptr = NULL;
		errno = 0;
		double const cut = strtod(s, &endptr);
		if (errno) {
			free(string);
			fprintf(stderr, "input: %s\n", strerror(errno));
			exit(EXIT_FAILURE);
		}

		if ((0 == cut) && (s == endptr)) {
			free(string);
			fprintf(stderr,
				"input: conversion error missing param ljcut\n");
			exit(EXIT_FAILURE);
		}

		errno = 0;
		s = endptr;
		endptr = NULL;
		double const offset = strtod(s, &endptr);
		if (errno) {
			free(string);
			fprintf(stderr, "input: %s\n", strerror(errno));
			exit(EXIT_FAILURE);
		}

		if ((0 == offset) && s == endptr) {
			free(string);
			fprintf(stderr,
				"input: conversion error missing param offsetflag\n");
			exit(EXIT_FAILURE);
		}

		if ((0 >= cut) || (0 != offset && 1 != offset)) {
			free(string);
			fprintf(stderr, "input: bad nonbond style param\n");
			exit(EXIT_FAILURE);
		}

		cutlj = cut;
		offsetflag = (int) offset;
		if (0 == mixflag) {
			mixstyle = 1;
		}

	} else {
		free(string);
		fprintf(stderr, "input: unsupported nonbond style\n");
		exit(EXIT_FAILURE);
	}

	if (0 == nonstyle) {
		fprintf(stdout, "input: Nonbond style none\n");
	} else if (1 == nonstyle) {
		fprintf(stdout,
			"input: Nonbond style %s %le %d\n",
			ljcutoff,
			cutlj,
			offsetflag);
	} else {
		free(string);
		fprintf(stderr, "input: implementation error\n");
		exit(EXIT_FAILURE);
	}
}

void input (int *iopflag)
{
	*iopflag = 0;
	int nlines = 0;
	while (1) {
		errno = 0;
		ssize_t rc = getline(&str, &sz_str, stdin);
		if (-1 == rc) {
			if (errno) {
				if (str) {
					free(str);
				}
				fprintf(stderr, "input: %s\n", strerror(errno));
				exit(EXIT_FAILURE);
			}
			*iopflag = EOF;
			break;
		}
		if (VERBOSE) {
			fprintf(stdout, "%s", str);
		}

		// units (if present) must be the first command in the input script

		if (strstr(str, "units")) {
			if (firstflag) {
				fprintf(stderr,
					"input: Units command must be the first command "
					"in file\n");
				exit(EXIT_FAILURE);
			}
			in_units(str);
			++nlines;
			continue;
		}

		firstflag = 1;

		// commands that must occur BEFORE read_data and read_restart

		if (strstr(str, "dimension")) {
			in_dimension(str);
			++nlines;
			continue;
		}

		// commands that can appear BEFORE or AFTER read_data and read_restart

		if (strstr(str, "nonbond style")) {
			in_nonbond_style(str);
			++nlines;
			continue;
		}

		// integrator settings

		if (strstr(str, "timestep")) {
			in_timestep(str);
			++nlines;
			continue;
		}
		++nlines;
	}
	if (VERBOSE) {
		fprintf(stdout, "nlines: %d\n", nlines);
	}
	if (str) {
		free(str);
	}
	sz_str = 0;
}
