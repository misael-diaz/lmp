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
// all global variables
// -------------------------------------------------------------------------
// *** parameters unlikely to need changing

#include <stddef.h>

#define MAXFIX 100
#define MAXDIAG 10

// maxfix = max # of fixes
// maxdiag = max # of user defined diagnostics

// minown = min size of atom arrays on one proc
// minghost = min size of ghost atom arrays on one proc
// minneigh = min size of neighbor list on one proc
// minbuf = min size of communication bufs on one proc

int const maxfix = MAXFIX;
int const maxdiag = MAXDIAG;

int const minown = 1000;
int const minghost = 1000;
int const minneigh = 10000;
int const minbuf = 10000;

// -------------------------------------------------------------------------
// *** constants

// boltz = Boltzmann factor
// dtfactor = multiplier on dt into time units
// efactor = conversion factor for Coulombic to energy units
// pfactor = conversion factor for pressure units
// two_1_3 = 2^(1/3)

double boltz = 0;
double dtfactor = 0;
double efactor = 0;
double pfactor = 0;
double two_1_3 = 0;

// -------------------------------------------------------------------------
// *** global settings and flags

// iversion = internal version # of code
// units = 0 for conventional, 1 for LJ units
// idimension = 2 or 3 for 2d/3d
// nsteps = # of timesteps to run
// itime = current timestep from 1 to nsteps
// ntimestep = current global timestep #
// ntime_last = what the global timestep will be on last step of run
// firstflag = 0 if haven't read any input command yet, 1 if have
// readflag = 0 if haven't read-in atoms yet, 1 if have
// newton = user input combination (0-3) of nonbond and bonded Newton flag
// newton_nonbond = 0 if Newton's law off for nonbond forces, 1 if on
// newton_bond = 0 if Newton's law off for bonded forces, 1 if on

// dt = timestep
// dthalf = 1/2 the timestep

int iversion = 0;
int units = 0;
int idimension = 0;
int nsteps = 0;
int itime = 0;
int ntimestep = 0;
int ntime_last = 0;
int firstflag = 0;
int readflag = 0;
int newton = 0;
int newton_nonbond = 0;
int newton_bond = 0;

double dt = 0;
double dthalf = 0;

// -------------------------------------------------------------------------
// *** domain

// xprd,yprd,zprd = size of global simulation box
// xprd_half,yprd_half,zprd_half = 1/2 the box lengths
// box(2,3) = lower/upper boundary of 3 dims of global box
// border(2,3) = lower/upper boundary of 3 dims of my sub-domain
// perflagx,perflagy,perflagz = 0 for periodic, 1 for non-periodic in 3 dims
// slabflag = 0 for 3-D periodicity of long-range forces, 1 for slab
// nonperiodic = TRUE if any dim is non-periodic, FALSE otherwise
// slab_volfactor = ratio of total volume to occupied volume for slab geometry
// zprd_slab = size of global box for slab geometry long-range interactions

double xprd = 0;
double yprd = 0;
double zprd = 0;
double xprd_half = 0;
double yprd_half = 0;
double zprd_half = 0;

double box[2][3];
double border[2][3];

int perflagx = 0;
int perflagy = 0;
int perflagz = 0;
int slabflag = 0;
int nonperiodic = 0;

double slab_volfactor = 0;
double zprd_slab = 0;

// -------------------------------------------------------------------------
// *** atoms

// natoms = total # of atoms in simulations
// ntypes = # of atom types
// maxown = most # of atoms I can own
// maxghost = most # of ghost atoms I can acquire
// extra_own = multiplier on allocation of maxown
// extra_ghost = multiplier on allocation of maxghost
// maxatom = maxown + maxghost
// maxatomp1 = maxatom + 1, used for storing a flag for special interactions

// nlocal = # of atoms I own
// nghost = # of ghost atoms I have
// max_nlocal = most atoms I ever own during simulation
// max_ghost = most ghost atoms I ever acquire during simulation

// mass(ntypes) = mass of each type of atom
// x(3,n) = coords of my atoms, owned and ghost
// v(3,n) = velocities of owned atoms
// f(3,n) = forces on my atoms, owned and ghost
// q(n) = charge of my atoms, owned and ghost
// xhold(3,n) = coords of owned atoms at last neighbor list build
// fhold(3,n) = force storage in integrate.f for newton=3 virial computation
// tag(n) = global ID tags of my atoms, owned and ghost
// type(n) = type of my atoms, owned and ghost
// molecule(n) = molecule ID # of owned atoms
// true(n) = true flag of owned atoms, which 3-d image of box they are in

int natoms = 0;
int ntypes = 0;
int maxown = 0;
int maxghost = 0;
double extra_own = 0;
double extra_ghost = 0;
int maxatom = 0;
int maxatomp1 = 0;

int nlocal = 0;
int nghost = 0;
int max_nlocal = 0;
int max_nghost = 0;

double *mass = NULL;

double *x = NULL;
double *v = NULL;
double *f = NULL;

double *q = NULL;
double *xhold = NULL;
double *fhold = NULL;
int *tag = NULL;
int *type = NULL;
int *molecule = NULL;
int *itrue = NULL;

// -------------------------------------------------------------------------
// *** bond connectivity for each atom

// nmolecular = 0 if an atomic system, 1 if molecular
// maxbondper,maxangleper,maxdihedper,maximproper = max # of bonds, etc
//       that any atom in simulation must store, depends on newton_bond
// maxspec = max # of 1-2, 1-3, 1-4 neighbors any atom must store

// numbond(n) = # of bonds stored by each of my atoms
// numangle(n) = # of angles stored by each of my atoms
// numdihed(n) = # of dihedrals stored by each of my atoms
// numimpro(n) = # of impropers stored by each of my atoms

// bondtype,angletype,dihedtype,improtype = 
//  type of each bond,angle,dihedral,improper for each of my atoms
// bondatom,angleatom,dihedatom,improatom 1234 = 
//  global tags of 2/3/4 atoms in each bond,angle,dihedral,improper
//  for each of my atoms

// num 123 bond = how many 1-2, 1-3, 1-4 neighbors this atom has in specbond
//                is cummulative: num2 includes 1-2 count, num3 has 1-2 & 1-3
// specbond = list of global ID tags of 
//            1-2, 1-3, 1-4 neighbors of each of my owned atoms

int nmolecular = 0;
int maxbondper = 0;
int maxangleper = 0;
int maxdihedper = 0;
int maximproper = 0;
int maxspec = 0;

int *numbond = NULL;
int *numangle = NULL;
int *numdihed = NULL;
int *numimpro = NULL;

int *bondtype = NULL;
int *bondatom1 = NULL;
int *bondatom2 = NULL;
int *angletype = NULL;
int *angleatom1 = NULL;
int *angleatom2 = NULL;
int *angleatom3 = NULL;
int *dihedtype = NULL;
int *dihedatom1 = NULL;
int *dihedatom2 = NULL;
int *dihedatom3 = NULL;
int *dihedatom4 = NULL;
int *improtype = NULL;
int *improatom1 = NULL;
int *improatom2 = NULL;
int *improatom3 = NULL;
int *improatom4 = NULL;

int *num1bond = NULL;
int *num2bond = NULL;
int *num3bond = NULL;
int *specbond = NULL;

// -------------------------------------------------------------------------
// *** global ptr

// localptr(i) = local index of atom with global tag of i
//               0 if this proc doesn't own the atom or have it as a ghost

int * localptr = NULL;

// -------------------------------------------------------------------------
// *** nonbond LJ and Coulombic interactions

// nonstyle = style of nonbond VanderWaal interactions
//            0 = none, 1 = cutoff LJ, 2 = smoothed LJ, 3 = shifted LJ,
//            4 = soft potential, 5 = class 2
// offsetflag = whether to add in shifted LJ energy at cutoff distance
// mixflag = 0 if no user input of mixing style, 1 if has been specified
// mixstyle = nonbond mixing style for epsilon,sigma
//            1 = geometric, 2 = arithmetic, 3 = sixthpower (class 2)
// ncharge = 0 if no charges defined in system, 1 if there are charges
// coulstyle = style of Coulombic interactions
//           0 = none, 1 = cutoff, 2 = smoothed, 3 = Ewald, 4 = PPPM
// amberflag = 0/1 if special bonds are set to AMBER force field settings
// coulpre = prefactor on Coulombic energy (efactor/dielectric)
// dielectric = dielectric constant settable by user
// special(3) = nonbond weighting factors on 1-2, 1-3, 1-4 neighbors

// cutforce,cutforce_sq = longest force cutoff of any nonbond interaction
// cutlj,cutlj_sq = cutoff for LJ
// cutljinterior,cutljint_sq = interior cutoff for smoothed LJ
// cutcoul = cutoff for Coulombic
// cutcoulsq = square of Coulombic cutoff
// cutcoulint = interior cutoff for smoothed Coulombic
// cutcoulintsq = square of interior cutoff for smoothed Coulombic
// cutmax = cutoff set by "maximum cutoff" command for memory allocator
// cutmemory = cutoff used by memory allocator, includes neighbor skin
// ch_denom_lj = charmm switching function denominator for LJ interactions
// ch_denom_coul = charmm switching fnct. denom. for Coulombic interactions
// kappa = damping factor for Debye/Huckel interactions

// noncoeff 1234 = nonbond coefficients as input for each atom type pair
// noncoeff14 12 = 1-4 interaction nonbond coefficients as input
// nontypeflag = whether coeffs have been specified for each type pair
// cutforcesq = force cutoff for each type pair
// cutneighsq = neighbor cutoff for each type pair
// cutljsq = LJ cutoff for each type pair
// cutljinner, cutljinnersq = interior LJ cutoff for each type pair
// lj 12345 = nonbond coeffs as used in force routines for each type pair
// lj14 1234 = 1-4 interaction nonbond coeffs as used in charmm force routine
// ljsw 01234 = smoothed LJ coeffs as used in force routines for each type pair
// offset = offset LJ energy for each type pair

int nonstyle = 0;
int offsetflag = 0;
int mixflag = 0;
int mixstyle = 0;
int ncharge = 0;
int coulstyle = 0;
int amberflag = 0;
double coulpre = 0;
double dielectric = 0;
double special[3];

double cutforce = 0;
double cutforce_sq = 0;
double cutlj = 0;
double cutljinterior = 0;
double cutlj_sq = 0;
double cutljint_sq = 0;
double cutcoul = 0;
double cutcoulsq = 0;
double cutcoulint = 0;
double cutcoulintsq = 0;
double cutmax = 0;
double cutmemory = 0;
double ch_denom_lj = 0;
double ch_denom_coul = 0;
double kappa = 0;

double *noncoeff1 = NULL;
double *noncoeff2 = NULL;
double *noncoeff3 = NULL;
double *noncoeff4 = NULL;
double *noncoeff14_1 = NULL;
double *noncoeff14_2 = NULL;
int *nontypeflag = NULL;

double *cutforcesq = NULL;
double *cutneighsq = NULL;
double *cutljsq = NULL;
double *cutljinner = NULL;
double *cutljinnersq = NULL;
double *lj1 = NULL;
double *lj2 = NULL;
double *lj3 = NULL;
double *lj4 = NULL;
double *lj5 = NULL;
double *lj14_1 = NULL;
double *lj14_2 = NULL;
double *lj14_3 = NULL;
double *lj14_4 = NULL;
double *ljsw0 = NULL;
double *ljsw1 = NULL;
double *ljsw2 = NULL;
double *ljsw3 = NULL;
double *ljsw4 = NULL;
double *offset = NULL;

// -------------------------------------------------------------------------
// *** bonded interactions

// nbonds = # of bonds in entire simulation
// nangles = # of angles in entire sim
// ndihedrals = # of dihedrals in entire sim
// nimpropers = # of impropers in entier sim
// nbondtypes,nangletypes,ndihedtypes,nimprotypes = 
//              # of types of each interaction
// bondstyle = style of bond interactions
//             0 = none, 1 = harmonic, 2 = FENE, 3 = shifted FENE,
//             4 = nonlinear spring, 5 = class 2
// anglestyle = style of angle interactions
//              0 = none, 1 = harmonic, 2 = class 2
// dihedstyle = style of dihedral interactions
//              0 = none, 1 = harmonic, 2 = class 2
// improstyle = style of improper interactions
//              0 = none, 1 = harmonic, 2 = cvff, 3 = class 2

// bondcoeff = coeffs for all bond types
// bondtypeflag = whether coeffs have been set for each type
// anglecoeff = coeffs for all angle types
// angletypeflag = whether coeffs have been set
// dihedcoeff = coeffs for all dihedral types
// dihedtypeflag = whether coeffs have been set
// improcoeff = coeffs for all improper types
// improtypeflag = whether coeffs have been set

int nbonds = 0;
int nangles = 0;
int ndihedrals = 0;
int nimpropers = 0;
int nbondtypes = 0;
int nangletypes = 0;
int ndihedtypes = 0;
int nimprotypes = 0;
int bondstyle = 0;
int anglestyle = 0;
int dihedstyle = 0;
int improstyle = 0;

double *bondcoeff = NULL;
int *bondtypeflag = NULL;
double *anglecoeff = NULL;
int *angletypeflag = NULL;
double *dihedcoeff = NULL;
int *dihedtypeflag = NULL;
double *improcoeff = NULL;
int *improtypeflag = NULL;

// -------------------------------------------------------------------------
// *** class 2 force field

// coeffs for all class 2 interactions

double *bondbondcoeff = NULL;
double *bondanglecoeff = NULL;
double *bondbond13coeff = NULL;
double *angleanglecoeff = NULL;
double *angleangletorsioncoeff = NULL;
double *angletorsioncoeff = NULL;
double *midbondtorsioncoeff = NULL;
double *endbondtorsioncoeff = NULL;

// -------------------------------------------------------------------------
// *** nonbond neighbor lists

// maxneigh = most # of neighbors that can be stored
// max_neigh = = most # ever stored during simulation
// extra_neigh = multiplier on allocation of maxneigh
// numneigh = # of times neighor lists are rebuilt
// ndanger = # of times a "dangerous" rebuild is done, where an atom may
//           have moved within a force-cutoff distance earlier
// maxbin = # of neighbor bins I own
// nbinx,nbiny,nbinz = # of bins in each dim of global simulation box
// mbinx,mbiny,mbinz = # of bins in each dim of my sub-domain, including ghosts
// mbinxlo,mbinylo,mbinzlo = lowest global bin one of my atoms could be in
// nstencil = # of bins in stencil for checking neighbor interactions
// stencil = list of offsets into set of bins that comprise the stencil
// binsizex,binsizey,binsizez = size of a neighbor bin in 3-d
// bininvx,bininvy,bininvz = inverse sizes of neighbor bins
// cutneigh = neighbor cutoff
// skin = distance beyond force cutoff that neighbor cutoff includes
// triggersq = trigger distance (1/2 of skin) for rebuilding neighbor lists
// neighago = how many steps ago neighbor list was rebuilt
// neighdelay = delay for at least this may steps before rebuilding lists
// neighfreq = check rebuild criterion or rebuild every this many steps
// neighstyle = 0 for brute-force N^2 search, 1 for binning
// neightrigger = 0 if don't rebuild based on distance moved, 1 if do
// nlist = 1-d list of all neighors of all my atoms
// nnfirst,nnlast = ptrs into nlist where neighbor of each atom start/stop
// bin = ptr from each atom to the next atom in its neighbor bin
// binpnt = ptr to 1st atom in each neighbor bin

int maxneigh = 0;
int max_neigh = 0;
double extra_neigh = 0;
int numneigh = 0;
int ndanger = 0;

int maxbin = 0;
int nbinx = 0;
int nbiny = 0;
int nbinz = 0;
int mbinx = 0;
int mbiny = 0;
int mbinz = 0;
int mbinxlo = 0;
int mbinylo = 0;
int mbinzlo = 0;

int nstencil = 0;
int stencil[1000];

double binsizex = 0;
double binsizey = 0;
double binsizez = 0;
double bininvx = 0;
double bininvy = 0;
double bininvz = 0;

double cutneigh = 0;
double skin = 0;
double triggersq = 0;

int neighago = 0;
int neighdelay = 0;
int neighfreq = 0;
int neighstyle = 0;
int neightrigger = 0;

int *nlist = NULL;
int *nnfirst = NULL;
int *nnlast = NULL;

int *bin = NULL;
int *binpnt = NULL;

// -------------------------------------------------------------------------
// *** bonded lists

// nbondlocal,nanglelocal,ndihedlocal,nimprolocal =
//   # of bond,angle,dihedral,improper in current bonded lists
// maxbondlocal,maxanglelocal,maxdihedlocal,maximprolocal
//   max # of bond,angle,dihedral,improper my allocated lists can store
// max_bond,max_angle,max_dihed,max_impro =
//   most # of bond,angle,dihedral,improper my lists every store
// bondlist,anglelist,dihedlist,improlist =
//   lists for bond,angle,dihedral,improper interactions to compute,
//   each entry in list stores type and global tag IDs of atoms involved

int nbondlocal = 0;
int nanglelocal = 0;
int ndihedlocal = 0;
int nimprolocal = 0;
int maxbondlocal = 0;
int maxanglelocal = 0;
int maxdihedlocal = 0;
int maximprolocal = 0;
int max_bond = 0;
int max_angle = 0;
int max_dihed = 0;
int max_impro = 0;

int *bondlist = NULL;
int *anglelist = NULL;
int *dihedlist = NULL;
int *improlist = NULL;

// -------------------------------------------------------------------------
// *** communication

// node = proc ID of me
// nprocs = total # of procs
// nswap = # of atom swaps each proc does with surrounding procs
// maxswap = current size of swap arrays
// max_exch = most # of atoms I've ever exchanged
// max_bord = most # of atoms sent in one border swap
// max_slist = most # of atoms I send in all my border swaps
//             (max_ghost is most I ever receive)
// maxbuf = size of allocated communication buffers
// max_buf = largest portion of buf ever used during simulation
// extra_buf = multiplier on allocation of maxbuf

// pgrid(3) = # of procs assigned to each dim on simulation box
// me(3) = which sub-box I own in each of 3 dims (0 to n-1)
// mpart(2,3) = my 6 neighboring procs
// need(3) = how many boxes away I need ghost info from in each dim

// slablo,slabhi = boundary inside which I need ghost info for in each swap
// spart = proc to send to in each swap
// rpart = proc to recv from in each swap
// nsfirst,nslast = ptrs into slist of send atoms for each swap
// nrfirst,nrlast = ptrs into ghost list where to put recv info in each swap
// commflag(3,n) = flags for PBC treatment of each dim in each swap
// commflagall(n) = flag for whether PBC treatment is needed in each swap
// slist = list of local atoms to send in all my swaps
// ibuf1,ibuf2,buf1,buf2 = buffers to use in communicating ghost info

int node = 0;
int nprocs = 0;
int nswap = 0;
int maxswap = 0;
int max_exch = 0;
int max_bord = 0;
int max_slist = 0;
int maxbuf = 0;
int max_buf = 0;
double extra_buf = 0;

int pgrid[3];
int me[3];
int mpart[2][3];
int need[3];

double *slablo = NULL;
double *slabhi = NULL;
int *spart = NULL;
int *rpart = NULL;
int *nsfirst = NULL;
int *nslast = NULL;
int *nrfirst = NULL;
int *nrlast = NULL;
int *commflag = NULL;
int *commflagall = NULL;

int *slist = NULL;

int *ibuf1 = NULL;
int *ibuf2 = NULL;
double *buf1 = NULL;
double *buf2 = NULL;

// -------------------------------------------------------------------------
// *** thermodynamics

// t_current = current temperature
// e_total = potential energy
// p_total = scalar pressure
// p_current(3) = pressure trace
// virial(6) = diagonal and off-diagonal virial components
// virialhold(6) = temporary copy of virial for ghost-atom virial computation
// vir_long(6) = virial for long-range Ewald/PPPM forces

double t_current = 0;
double e_total = 0;
double p_total = 0;
double p_current[3];
double virial[6];
double virialhold[6];
double vir_long[6];

// nonbond and bonded potential energies
// e_14 terms are for nonbond portion computed in CHARMM dihedrals

double e_vdwl = 0;
double e_coul = 0;
double e_bond = 0;
double e_angle = 0;
double e_dihedral = 0;
double e_improper = 0;
double e_14_coul = 0;
double e_14_vdwl = 0;

// long-range Coulombic energy

double e_long = 0;

// class 2 potential energies

double e_bondbond = 0;
double e_bondbond13 = 0;
double e_bondangle = 0;
double e_endbondtorsion = 0;
double e_midbondtorsion = 0;
double e_angletorsion = 0;
double e_angleangletorsion = 0;
double e_angleangle = 0;

// -------------------------------------------------------------------------
// *** ensemble variables: temp, pressure, volume control

// ensemble = which ensemble is being simulated: 1=NVE, 2=NVT, 3=NPH, 4=NPT

// tempstyle = style of temperature control
//             0 = none, 1 = rescale, 2 = replace, 3 = Langevin, 4 = Nose/Hoover
// t_every = check for temp rescaling every this many steps
// t_start,t_stop = desired temperature at start/end of run
// t_window = rescale temperature if it is outside this window
// t_fraction = amount (0% to 100%) of rescaling to perform
// t_freq = drag/mass parameter in Langevin and Nose/Hoover temp control
// t_target = target temperature on this timestep
// eta,eta_dot = evolving Nose variable for temp control

// pressstyle = style of pressure control: 0 = none, 1 = Nose/Hoover
// presscouple = style of coupling:
//               0 = xyz isotropic, 1 = xy, 2 = yz, 3 = xz, 4 = anisotropic
// p_freq(3) = piston mass parameter in Nose/Hoover pressure control
// p_start(3),p_stop(3) = desired pressure at start/end of run
// p_target(3) = target pressure on this timestep
// omega(3),omega_dot(3) = evolving Nose variables for pressure control
// masssum = total mass in system

// volstyle = style of volume contfol: 0 = none, 1 = linear expand/contract
// voldimx,voldimy,voldimz = disabled (0) or enabled (1) for each dimension
// volstart/stop xyz lo/hi = global box boundary at begin/end of run in each dim

int ensemble = 0;

int tempstyle = 0;
int t_every = 0;
double t_start = 0;
double t_stop = 0;
double t_window = 0;
double t_fraction = 0;
double t_freq = 0;
double t_target = 0;
double eta = 0;
double eta_dot = 0;

int pressstyle = 0;
int presscouple = 0;
double p_freq[3];
double p_start[3];
double p_stop[3];
double p_target[3];
double omega[3];
double omega_dot[3];
double masssum = 0;

int volstyle = 0;
int voldimx = 0;
int voldimy = 0;
int voldimz = 0;
double volstart_xlo = 0;
double volstart_xhi = 0;
double volstop_xlo = 0;
double volstop_xhi = 0;
double volstart_ylo = 0;
double volstart_yhi = 0;
double volstop_ylo = 0;
double volstop_yhi = 0;
double volstart_zlo = 0;
double volstart_zhi = 0;
double volstop_zlo = 0;
double volstop_zhi = 0;

// -------------------------------------------------------------------------
// *** temperature creation

// createstyle = style of velocity creation
//               1 = uniform, 2 = gaussian, 3 = explicit velocity
// creategroup = kind of group of atoms to create vels for
//               1 = types, 2 = region, 3 = remainder of unset atoms
// createlo,createhi = range of types/molecules to create vels for
// iseed = random # seed, used for velocity creation and temp control
// rotationflag = 0 if don't zero angular momentum, 1 if do
// t_create = desired temperature to create vels for
// createregion(6) = geometric bounds of atoms to create vels for
// createvec(3) = explicit velocity vector to apply to atoms in group

// velflag(:) = used to mark atoms with previously initialized vels

int createstyle = 0;
int creategroup = 0;
int createlo = 0;
int createhi = 0;
int iseed = 0;
int rotationflag = 0;
double t_create = 0;
double createregion[6];
double createvec[3];

int* velflag = NULL;

// -------------------------------------------------------------------------
// *** fixes

// nfixes = # of fixes specified by user
// nfixes_respa = # of setforce and aveforce fixes that must be applied
//                at each level of rRESPA

// fixstyle = style of each fix, 0 = none, 1 = setforce, 2 = addforce, 
//            3 = aveforce, 4 = rescale, 5 = Langevin, 6 = Nose/Hoover,
//            7 = springforce, 8 = dragforce, 9 = shake
// fixflag = flags associated with each fix
// fixptr = ptr into fixstore where values associated with a fix are stored
// fixcount = # of atoms assigned to each fix
// fixactive = flag for whether each fix is active on this timestep
// fixcoeff = coeffs/values/constants assocated with each fix
// fixmass = total mass of all atoms assigned to each fix

// fixwhich = which fix (1-N) to assign atoms to
// fixgroup = kind of group that the fix is to be applied to
//            1 = single atom, 2 = molecule, 3 = type, 4 = region, 5 = remainder
// fixatom = atom tag ID or molecule ID to apply fix to
// fixtype,fixbond,fixangle = atom/bond/angle type to apply fix to
// fixregion(6) = geometric region to apply fix to

// fixnum = how many fix quantities need to be summed up across procs
// fixstore = vector of fix values computed on this timestep
// fixstoretmp,fixmasstmp,fixcounttmp = scratch vectors for summing fix info

// fix = fix #'s associated with each of my atoms

int nfixes = 0;
int nfixes_respa = 0;

int fixstyle[MAXFIX];
int fixflag[3][MAXFIX];
int fixptr[MAXFIX];
int fixcount[MAXFIX];
int fixactive[MAXFIX];
double fixcoeff[7][MAXFIX];
double fixmass[MAXFIX];

int fixwhich = 0;
int fixgroup = 0;
int fixatom = 0;
int fixtype = 0;
int fixbond = 0;
int fixangle = 0;
double fixregion[6];

int fixnum = 0;
double fixstore[3*MAXFIX];

double fixstoretmp[3*MAXFIX];
double fixmasstmp[MAXFIX];
int fixcounttmp[MAXFIX];

int *fix = NULL;

// -------------------------------------------------------------------------
// *** Ewald (some of these variables are also used by PPPM)

// kcount = actual # of Ewald vectors
// kmax = max dimensionality of Ewald vectors
// gewald = G vector for Ewald/PPPM as function of precision/cutoff
// qsum,qsqsum = total charge, square of total charge
// long_prec = user-specified precision for long-range approximation

// kxvecs,kyvecs,kzvecs = pre-computed indices of Ewald vectors
// ug = pre-computed exponential factor for each Ewald vector
// eg(3,n) = pre-computed Ewald vector components of ug
// vg(6,n) = pre-computed virial components for each Ewald/PPPM vector
// sfacrl,sfacim = real/imag structure factors for each Ewald vector
// sfacrl_all,sfacim_all = scratch vectors for summing sfacs across procs
// cs,sn = atom contributions to each Ewald level

// ek(3,n) = components of electric field at each atom

int kcount = 0;
int kmax = 0;
double gewald = 0;
double qsum = 0;
double qsqsum = 0;
double long_prec = 0;

int *kxvecs = NULL;
int *kyvecs = NULL;
int *kzvecs = NULL;
double *ug = NULL;
double *eg = NULL;
double *vg = NULL;
double *sfacrl = NULL;
double *sfacim = NULL;
double *sfacrl_all = NULL;
double *sfacim_all = NULL;
double *cs = NULL;
double *sn = NULL;

double *ek = NULL;

// -------------------------------------------------------------------------
// *** PPPM

// orderflag = order of PPPM, how far into grid the charge overlaps
// meshflag = 1 if user sets PPPM mesh, 0 otherwise
// nfft = # of FFT points I own in FFT decomp
// nlower,nupper = # of grid pts a charge extends to the left,right
// nx_pppm,ny_pppm,nz_pppm = global PPPM grid
// nx_pppm_input,ny_pppm_input,nz_pppm_input = global PPPM grid chosen by user
// nxlo_in,nxhi_in,nylo_in,nyhi_in,nzlo_in,nzhi_in = 
//   portion of global grid I own in brick decomp
// nxlo_out,nxhi_out,nylo_out,nyhi_out,nzlo_out,nzhi_out = 
//   portion of global grid I own including ghost cells
// nxlo_ghost,nxhi_ghost,nylo_ghost,nyhi_ghost,nzlo_ghost,nzhi_ghost =
//   # of planes of ghost cells I receive from neighbor in each direction
// nxlo_fft,nxhi_fft,nylo_fft,nyhi_fft,nzlo_fft,nzhi_fft =
//   portion of global grid I own in FFT decomp
// plan1_fft,plan2_fft,plan_remap = FFT and remap plans

// maxgrid = size of local grid array including ghost cells in brick decomp
// maxfft = size of local FFT array, bigger of brick or FFT decomp
// maxpbuf = size of buffers for exchanging ghost cells

// density_brick = stores particle charge as mapped to grid, brick decomp
// vdx_brick,vdy_brick,vdz_brick = stores potential gradient as mapped to grid
// density_fft = particle charge on grid in FFT decomp
// greensfn = Green's function for each point on grid
// workvec1,workvec2 = complex work vectors for FFTs
// partgrid = which grid cell (nx,ny,nz) a particle is centered at
// fkvecs_x,fkvecs_y,fkvecs_z = pre-computed FFT coefficients
// pbuf1,pbuf2 = buffers for exchanging ghost cells

int orderflag = 0;
int meshflag = 0;
int nfft = 0;
int nlower = 0;
int nupper = 0;
int nx_pppm = 0;
int ny_pppm = 0;
int nz_pppm = 0;
int nx_pppm_input = 0;
int ny_pppm_input = 0;
int nz_pppm_input = 0;
int nxlo_in = 0;
int nxhi_in = 0;
int nylo_in = 0;
int nyhi_in = 0;
int nzlo_in = 0;
int nzhi_in = 0;
int nxlo_out = 0;
int nxhi_out = 0;
int nylo_out = 0;
int nyhi_out = 0;
int nzlo_out = 0;
int nzhi_out = 0;
int nxlo_ghost = 0;
int nxhi_ghost = 0;
int nylo_ghost = 0;
int nyhi_ghost = 0;
int nzlo_ghost = 0;
int nzhi_ghost = 0;
int nxlo_fft = 0;
int nxhi_fft = 0;
int nylo_fft = 0;
int nyhi_fft = 0;
int nzlo_fft = 0;
int nzhi_fft = 0;
double plan1_fft = 0;
double plan2_fft = 0;
double plan_remap = 0;

int maxgrid = 0;
int maxfft = 0;
int maxpbuf = 0;

double *density_brick = NULL;
double *vdx_brick = NULL;
double *vdy_brick = NULL;
double *vdz_brick = NULL;
double *density_fft = NULL;
double *greensfn = NULL;
/*
complex*16 *workvec1(:),workvec2(:)
*/
int *partgrid = NULL;
double *fkvecs_x = NULL;
double *fkvecs_y = NULL;
double *fkvecs_z = NULL;
double *pbuf1 = NULL;
double *pbuf2 = NULL;

// -------------------------------------------------------------------------
// *** rRESPA

// nrespa = 1 if rRESPA is on, 0 if off
// nstretch,nintra,nshort = dilation factors between 4 sets of timesteps
// vir_stretch(6),vir_intra(6),vir_short(6) = virials for 3 forces
// dthalf_intra,dthalf_short,dthalf_long = timesteps for 3 scales besides dt
// f_stretch,f_intra,f_short,f_long = force arrays for 4 forces

int nrespa = 0;
int nstretch = 0;
int nintra = 0;
int nshort = 0;
double vir_stretch[6];
double vir_intra[6];
double vir_short[6];
double dthalf_intra = 0;
double dthalf_short = 0;
double dthalf_long = 0;

double *f_stretch = NULL;
double *f_intra = NULL;
double *f_short = NULL;
double *f_long = NULL;

// -------------------------------------------------------------------------
// *** SHAKE

// nshake = 1 if SHAKE is on, 0 if off
// shakeiter = max # of iterations SHAKE will attempt
// nshakestats = print bond statistics every this many timesteps (0 = never)
// shakewhichbondcoeff = which bond coeff is bond length (for bondstyle)
// shakeableangle = 0/1 if angle constraint is not set or set
// shakeableanglebond = which bond type the angle constraint applies to
// nshake_next = next timestep to call bond statistics routine
// shaketol = tolerance for SHAKE bonds
// shakeanglebond = psuedo-bond distance across constrained angle
// shakeablebond = 0/1 for each bondtype, whether (0) not SHAKE or (1) SHAKE
// nshakebonds = total # of bonds & pseudo-bonds constrained by SHAKE

// shakegroup = 0,2,3,4 = size of SHAKE group this atom is part of
// shakepartner = global ID tag of all atoms (including self) in SHAKE group,
//                central atom is listed 1st
// shakebondtype = which type of bond each SHAKE bond is

// nshakelocal = # of SHAKE groups this proc must currently compute (neigh list)
// shakesize = size of each SHAKE group (2,3,4) (in neigh list)
// shakeatom = local IDs for each atom in SHAKE group (in neigh list)
// shakebondlen = bond length of each bond in SHAKE group (in neigh list)
// xshake = updated unconstrained atom coords for owned and ghost atoms

int nshake = 0;
int shakeiter = 0;
int nshakestats = 0;
int nshake_next = 0;
int nshakebonds = 0;
int shakewhichbondcoeff = 0;
int shakeableangle = 0;
int shakeableanglebond = 0;
double shaketol = 0;
double shakeanglebond = 0;
int *shakeablebond = NULL;

int *shakegroup = NULL;
int *shakepartner = NULL;
int *shakebondtype = NULL;

int nshakelocal = 0;
int *shakesize = NULL;
int *shakeatom = NULL;
double *shakebondlen = NULL;
double *xshake = NULL;

// -------------------------------------------------------------------------
// *** minimizer

// optstyle = style of minimizer, 1 = Hessian-free truncated Newton
// optflag = output minimizer iteration info every this many steps
// opt_max_iters,opt_max_fns = max # of iterations and function evals for min
// opt_stop_tol = stopping tolerance

int optstyle = 0;
int optflag = 0;
int opt_max_iters = 0;
int opt_max_fns = 0;
double opt_stop_tol = 0; 

// -------------------------------------------------------------------------
// *** output

// noutput_next = next timestep to do any kind of output on
// trueflag = flag for whether true box flags are included on input/output

// thermostyle = style of thermo output, little to lots
// nthermo,nthermo_next = how-often,when-next to do thermo output
// ndumpatom,ndumpatom_prev,ndumpatom_next = how-often,when-next to dump atoms
// ndumpvel,ndumpvel_prev,ndumpvel_next = how-often,when-next to dump vels
// ndumpforce,ndumpforce_prev,ndumpforce_next = how-often,when-next dump forces
// dumpatomfileflag,dumpvelfileflag,dumpforcefileflag = flags for whether
//   dump files are currently open or closed
// nrestart,nrestart_next = how-often,when-next to write a restart file
// restartstyle = (1) append timestep to file, (2) toggle between 2 files
// restartlast = which file (1 or 2) was used on last restart output
// restart_version = what version to expect when reading in a restart file

int noutput_next = 0;
int trueflag = 0;

int thermostyle = 0;
int nthermo = 0;
int nthermo_next = 0;
int ndumpatom = 0;
int ndumpatom_prev = 0;
int ndumpatom_next = 0;
int ndumpvel = 0;
int ndumpvel_prev = 0;
int ndumpvel_next = 0;
int ndumpforce = 0;
int ndumpforce_prev = 0;
int ndumpforce_next = 0;
int dumpatomfileflag = 0;
int dumpvelfileflag = 0;
int dumpforcefileflag = 0;
int nrestart = 0;
int nrestart_next = 0;
int restartstyle = 0;
int restartlast = 0;
int restart_version = 0;

// files for various I/O operations

char datafile[80];
char dumpatomfile[80];
char dumpvelfile[80];
char dumpforcefile[80];
char restart_in[80];
char restart_out[80];
char restart_out1[80];
char restart_out2[80];

// -------------------------------------------------------------------------
// *** diagnostics

// numdiag = # of diagnostic routines defined by user in diagnostic.f
// diagnames = name of each diag routine as specified in input script

// ndiag_next = next timestep to call any diagnostic routine on
// ndiag = how often to call each diag routine
// diagprev,diagnext = previous/next timestep to call each diag routine
// diagcall = flag for whether each diag routine should be called this timestep

// diagnparams = # of params to pass to each diag routine
// diagparam(5,n) = parameters to pass to each diag routine
// diagfileflag = whether file for each diag routine is open or closed
// diagfile = filename for each diag routine

int numdiag = 0;
char diagnames[MAXDIAG][16];

int ndiag_next = 0;
int ndiag[MAXDIAG];
int diagprev[MAXDIAG];
int diagnext[MAXDIAG];
int diagcall[MAXDIAG];

int diagnparams[MAXDIAG];
double diagparam[MAXDIAG][5];
int diagfileflag[MAXDIAG];
char diagfile[MAXDIAG][80];

// -------------------------------------------------------------------------
// *** timers

// CPU timers for various operations within code

double time_total = 0;
double time_loop = 0;
double time_current = 0;
double time_nonbond = 0;
double time_bond = 0;
double time_angle = 0;
double time_dihedral = 0;
double time_improper = 0;
double time_comm = 0;
double time_fcomm = 0;
double time_exch = 0;
double time_io = 0;
double time_shake = 0;
double time_neigh1 = 0;
double time_neigh2 = 0;
double time_long = 0;
double time_rho = 0;
double time_poiss = 0;
double time_field = 0;
double time_other = 0;

// -------------------------------------------------------------------------
