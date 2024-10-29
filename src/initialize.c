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
// initialize run parameters
// -------------------------------------------------------------------------

extern int iversion;

extern double boltz;
extern double dtfactor;
extern double pfactor;
extern double efactor;
extern double two_1_3;

extern int units;
extern int idimension;
extern int nonperiodic;
extern double skin;

extern double dt;

extern int volstyle;
extern int tempstyle;
extern int pressstyle;

extern int nonstyle;
extern double cutlj;
extern int coulstyle;

extern int nfixes;

extern int firstflag;
extern int readflag;

void initialize (void)
{
	iversion = 2001;

	// constants - init default is conventional units

	boltz = 0.001987191;
	dtfactor = 48.88821;
	pfactor = 68589.796;
	efactor = 332.0636;
	two_1_3 = 1.2599210498948732;

	// values that must be set before read_data or read_restart commands

	units = 0;
	idimension = 3;
	nonperiodic = 0;
	skin = 2.0;

	// integrator settings

	dt = 1.0;

	// ensemble settings

	volstyle = 0;
	tempstyle = 0;
	pressstyle = 0;

	// forcefield settings

	nonstyle = 1;
	cutlj = 10.0;
	coulstyle = 1;

	// no fixes

	nfixes = 0;

	// no commands have been read-in yet

	firstflag = 0;
	readflag = 0;
}
