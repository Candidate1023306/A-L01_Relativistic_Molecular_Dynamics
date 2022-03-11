// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

//DEC 20 2021, added gamma and momentum to the atom.cpp and atom.h files in the current build. Let's see if they fuck up or not. 



#include "fix_nve.h"
#include <iostream>

#include "atom.h"
#include "error.h"
#include "force.h"
#include "respa.h"
#include "update.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixNVE::FixNVE(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (!utils::strmatch(style,"^nve/sphere") && narg < 3)
    error->all(FLERR,"Illegal fix nve command");

  dynamic_group_allow = 1;
  time_integrate = 1;
}

/* ---------------------------------------------------------------------- */

int FixNVE::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  mask |= INITIAL_INTEGRATE_RESPA;
  mask |= FINAL_INTEGRATE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixNVE::init()
{
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;

  if (utils::strmatch(update->integrate_style,"^respa"))
    step_respa = ((Respa *) update->integrate)->step;
}

/* ----------------------------------------------------------------------
   allow for both per-type and per-atom mass
------------------------------------------------------------------------- */

void FixNVE::initial_integrate(int /*vflag*/)
{
  double dtfm;

  // update v and x of atoms in group

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double c = 2.998e+6;
  double c_squared = pow(c, 2);
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  if (rmass) {
    for (int i = 0; i < nlocal; i++) //i is the atom index. 
      if (mask[i] & groupbit) {
        //this is the bit I am going to start messing with. Commenting out momentum version for now jan 8 22
        dtfm = dtf / 39.95; //argon mass. 
        double v_squared = v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2];
        double gamma = 1/sqrt(1-v_squared/c_squared);
        // implementing relativistic equations of motion
        // updating velocity
        double big_gamma_x = ((gamma*v[i][0]) + (dtfm*f[i][0]));
        double big_gamma_y = ((gamma*v[i][1]) + (dtfm*f[i][1]));
        double big_gamma_z = ((gamma*v[i][2]) + (dtfm*f[i][2]));
        double big_gamma_dot = big_gamma_x*big_gamma_x + big_gamma_y*big_gamma_y + big_gamma_z*big_gamma_z; //get rid of pow. 
        //swap over sqrt to a new line.   
        v[i][0] = big_gamma_x/sqrt(1 + big_gamma_dot/c_squared);
        v[i][1] = big_gamma_y/sqrt(1 + big_gamma_dot/c_squared);
        v[i][2] = big_gamma_z/sqrt(1 + big_gamma_dot/c_squared);
        // updating position
        x[i][0] += dtv * v[i][0];
        x[i][1] += dtv * v[i][1];
        x[i][2] += dtv * v[i][2];
            }

  } else {
    for (int i = 0; i < nlocal; i++) //don't care about these, as I am using rmass. 
      if (mask[i] & groupbit) {
        dtfm = dtf / mass[type[i]];
        v[i][0] += dtfm * f[i][0];
        v[i][1] += dtfm * f[i][1];
        v[i][2] += dtfm * f[i][2];
        x[i][0] += dtv * v[i][0];
        x[i][1] += dtv * v[i][1];
        x[i][2] += dtv * v[i][2];
      }
  }
}

/* ---------------------------------------------------------------------- */

void FixNVE::final_integrate() //forces have been updated in the meantime. 
{
  double dtfm;

  // update v of atoms in group

  double **v = atom->v;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double c = 2.998e+6;
  double c_squared = pow(c,2);
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        dtfm = dtf / 39.95;
        double v_squared = v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2];
        double gamma = 1/sqrt(1-v_squared/c_squared);
        // implementing relativistic equations of motion
        // updating velocity
        double big_gamma_x = ((gamma*v[i][0]) + (dtfm*f[i][0]));
        double big_gamma_y = ((gamma*v[i][1]) + (dtfm*f[i][1]));
        double big_gamma_z = ((gamma*v[i][2]) + (dtfm*f[i][2]));
        double big_gamma_dot = big_gamma_x*big_gamma_x + big_gamma_y*big_gamma_y + big_gamma_z*big_gamma_z;  
        v[i][0] = big_gamma_x/sqrt(1 + big_gamma_dot/c_squared);
        v[i][1] = big_gamma_y/sqrt(1 + big_gamma_dot/c_squared);
        v[i][2] = big_gamma_z/sqrt(1 + big_gamma_dot/c_squared);
            }

  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        dtfm = dtf / mass[type[i]];
        v[i][0] += dtfm * f[i][0];
        v[i][1] += dtfm * f[i][1];
        v[i][2] += dtfm * f[i][2];
      }
  }
}

/* ---------------------------------------------------------------------- */

void FixNVE::initial_integrate_respa(int vflag, int ilevel, int /*iloop*/)
{
  dtv = step_respa[ilevel];
  dtf = 0.5 * step_respa[ilevel] * force->ftm2v;

  // innermost level - NVE update of v and x
  // all other levels - NVE update of v

  if (ilevel == 0) initial_integrate(vflag);
  else final_integrate();
}

/* ---------------------------------------------------------------------- */

void FixNVE::final_integrate_respa(int ilevel, int /*iloop*/)
{
  dtf = 0.5 * step_respa[ilevel] * force->ftm2v;
  final_integrate();
}

/* ---------------------------------------------------------------------- */

void FixNVE::reset_dt()
{
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
}