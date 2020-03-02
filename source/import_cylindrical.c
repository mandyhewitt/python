/***********************************************************/
/** @file  import_cylindrical.c
 * @author ksl
 * @date   May, 2018
 *
 * @brief  Routines to read in an arbitrary wind model in
 * cylindrical coordinates.
 *
 * These routines have been tested with models for FU Ori
 * produced by Lee Hartmann
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"
#include "import.h"


/**********************************************************/
/**
 * @brief      Read the an arbitrary wind model in cylindrical
 *     coordinates
 *
 * @param [in] ndom   The domain number for the imported model
 * @param [in] filename   The file containing the model to import
 * @return     Always returns 0
 *
 * @details
 *
 * This routine just reads in the data and stores it in arrays
 *
 * ### Notes ###
 * The basic data we need to read in are
 *
 * * i, j, inwind, x z  v_x v_y v_z rho (and optionally T)
 *
 * where v_x,_v_y,v_z are the velocities in the x,z plane
 * and where
 *
 * * i is the column number  (Thus i corresponds to ndim)
 * * j is the row number     (and z coresponds to mdim)
 * * inwind indicates whether this cell is in the wind
 *
 * We assume that all of the variables are defined as they
 * are in the python structues, that is that positions and
 * velocities are defined at the corners of cells but that
 * rho refers to the cell center.
 *
 * Note that we assume that the data are being read in in the
 * same order as printed out by windsave2table, that is that
 * the first "column" is read in, and then the second "column".
 * The grid must also be complete. Missing i or j elements are
 * not allowed.
 *
 **********************************************************/

int
import_cylindrical (ndom, filename)
     int ndom;
     char *filename;
{
  FILE *fptr;
  char line[LINELENGTH];
  int n, icell, jcell, ncell, inwind;
  double x, z, v_x, v_y, v_z, rho, t_r, t_e;
  int jz, jx;
  double delta;

  Log ("Reading a model in cylindrical coordinates %s\n", filename);

  if ((fptr = fopen (filename, "r")) == NULL)
  {
    Error ("import_cylindrical: No such file\n");
    Exit (1);                   /* No need to worry about mp at this point */
  }

  ncell = 0;
  while (fgets (line, LINELENGTH, fptr) != NULL)
  {
    n = sscanf (line, " %d %d %d %le %le %le %le %le %le %le %le", &icell, &jcell, &inwind, &x, &z, &v_x, &v_y, &v_z, &rho, &t_r, &t_e);

    if (n < READ_NO_TEMP_2D)
    {
      continue;
    }
    else
    {
      imported_model.i[ncell] = icell;
      imported_model.j[ncell] = jcell;
      imported_model.inwind[ncell] = inwind;
      imported_model.x[ncell] = x;
      imported_model.z[ncell] = z;
      imported_model.v_x[ncell] = v_x;
      imported_model.v_y[ncell] = v_y;
      imported_model.v_z[ncell] = v_z;
      imported_model.mass_rho[ncell] = rho;

      if (n >= READ_RAD_TEMP_2D)
      {
        imported_model.t_r[ncell] = t_r;
      }
      else if (n >= READ_BOTH_TEMP_2D)
      {
        imported_model.t_r[ncell] = t_r;
        imported_model.t_e[ncell] = t_e;
      }
      else
      {
        imported_model.t_r[ncell] = DEFAULT_IMPORT_TEMPERATURE;
      }

      ncell++;

      if (ncell > NDIM_MAX * NDIM_MAX)
      {
        Error ("%s : %i : trying to read in more grid points than allowed (%i). Try changing NDIM_MAX and recompiling.\n", __FILE__,
               __LINE__, NDIM_MAX * NDIM_MAX);
        Exit (1);
      }

    }
  }


/* 
 * Now calculate the dimensions of the grid.  This next calculation makes the assumption that
 * The last element of the grid was the last grid cell.  So we now calculate the sizes of
 * the grid.
 *
 * Although the initialization of most of zdom should be postponed
 * one has to give zdom the dimensions of the array; otherwise
 * the wrong number of elements in wmain for the wind will be allocated
 */

  zdom[ndom].ndim = imported_model.ndim = icell + 1;
  zdom[ndom].mdim = imported_model.mdim = jcell + 1;
  zdom[ndom].ndim2 = zdom[ndom].ndim * zdom[ndom].mdim;
  imported_model.ncell = ncell;

  /* Check that the grid is complete */

  if (ncell != imported_model.ndim * imported_model.mdim)
  {
    Error ("The dimensions of the imported grid seem wrong % d x %d != %d\n", imported_model.ndim, imported_model.mdim,
           imported_model.ncell);
    Exit (1);
  }

/* Fill 1-d arrays in x and z */

  jz = jx = 0;
  for (n = 0; n < imported_model.ncell; n++)
  {
    if (imported_model.i[n] == 0)
    {
      imported_model.wind_z[jz] = imported_model.z[n];
      jz++;
    }
    if (imported_model.j[n] == 0)
    {
      imported_model.wind_x[jx] = imported_model.x[n];
      jx++;
    }
  }


  /* Now fill in wind_midx and midz. Given how we construct
   * the mdpts we need to add one more on the end, though it
   * is not entirely obvious that this is needed, given the
   * assumption that we do not need extra buffer cells */

  for (n = 0; n < jz - 1; n++)
  {
    imported_model.wind_midz[n] = 0.5 * (imported_model.wind_z[n] + imported_model.wind_z[n + 1]);
  }


  delta = (imported_model.wind_z[n - 1] - imported_model.wind_z[n - 2]);
  imported_model.wind_midz[n] = imported_model.wind_z[n - 1] + 0.5 * delta;



  for (n = 0; n < jx - 1; n++)
  {
    imported_model.wind_midx[n] = 0.5 * (imported_model.wind_x[n] + imported_model.wind_x[n + 1]);
  }



  delta = (imported_model.wind_x[n - 1] - imported_model.wind_x[n - 2]);
  imported_model.wind_midx[n] = imported_model.wind_x[n - 1] + 0.5 * delta;


  return (0);
}




/**********************************************************/
/**
 * @brief       Use the imported data to initialize various
 * portions of the Wind and Domain structures
 *
 *
 * @param [in] w   The entire wind
 * @param [in] ndom   The domain number
 * @return   Always returns 0
 *
 * @details
 * This routine initializes the portions of the wind structure
 * using the imported model, specifically those portions having
 * to do with positions, and velocities.
 *
 * The routine creates a pillbox around the grid to be used
 * for defining the region which the wind (maximally) occupies.
 *
 *
 * ### Notes ###
 *
 **********************************************************/

int
cylindrical_make_grid_import (w, ndom)
     WindPtr w;
     int ndom;
{
  int n;
  int nn;
  double rmin, rmax, rho_min, rho_max, zmin, zmax;
  double x[3], r_inner, r_outer;

/*  This is an attempt to make the grid directly.
 *
 *  Note also that none of this will work unless a complete grid is read
 *  in
 *  */
  for (n = 0; n < imported_model.ncell; n++)
  {
    wind_ij_to_n (ndom, imported_model.i[n], imported_model.j[n], &nn);
    w[nn].x[0] = imported_model.x[n];
    w[nn].x[1] = 0;
    w[nn].x[2] = imported_model.z[n];
    w[nn].v[0] = imported_model.v_x[n];
    w[nn].v[1] = imported_model.v_y[n];
    w[nn].v[2] = imported_model.v_z[n];
    w[nn].inwind = imported_model.inwind[n];

    if (w[nn].inwind == W_NOT_INWIND || w[nn].inwind == W_PART_INWIND)
      w[nn].inwind = W_IGNORE;

    w[nn].xcen[0] = imported_model.wind_midx[imported_model.i[n]];
    w[nn].xcen[1] = 0;
    w[nn].xcen[2] = imported_model.wind_midz[imported_model.j[n]];

  }



  /* Now add information used in zdom */

  for (n = 0; n < zdom[ndom].ndim; n++)
  {
    zdom[ndom].wind_x[n] = imported_model.wind_x[n];
  }



  for (n = 0; n < zdom[ndom].mdim; n++)
  {
    zdom[ndom].wind_z[n] = imported_model.wind_z[n];
  }


  /* Now set up wind boundaries so they are harmless.

   * Note that given that we have already filled out
   * the WindPtr it seems like we could do this
   * without needing the internal structure used in this 
   * routine.  It's done this way to follow the procedure
   * in the rtheta model
   * 
   * Note that one has to be careful here, because one
   * has to include all of the outer most cell that
   * is in the wind.
   * 
   */


  rmax = rho_max = zmax = 0;
  rmin = rho_min = zmin = VERY_BIG;
  for (n = 0; n < imported_model.ncell; n++)

  {
    if (imported_model.inwind[n] >= 0)
    {
      x[0] = imported_model.x[n];
      x[1] = 0;
      x[2] = imported_model.z[n];

      r_inner = length (x);

      x[0] = imported_model.x[n + imported_model.mdim];
      x[1] = 0;
      x[2] = imported_model.z[n + 1];

      r_outer = length (x);

      if (imported_model.x[n + imported_model.mdim] > rho_max)
      {
        rho_max = imported_model.x[n + imported_model.mdim];
      }
      if (imported_model.z[n + 1] > zmax)
      {
        zmax = imported_model.z[n + 1];
      }
      if (imported_model.z[n] < zmin)
      {
        zmin = imported_model.z[n];
      }
      if (r_outer > rmax)
      {
        rmax = r_outer;
      }
      if (rho_min > imported_model.x[n])
      {
        rho_min = imported_model.x[n];
      }
      if (rmin > r_inner)
      {
        rmin = r_inner;
      }
    }
  }



  Log ("Imported:    rmin    rmax  %e %e\n", rmin, rmax);
  Log ("Imported:    zmin    zmax  %e %e\n", zmin, zmax);
  Log ("Imported: rho_min rho_max  %e %e\n", rho_min, rho_max);

  zdom[ndom].wind_rho_min = zdom[ndom].rho_min = rho_min;
  zdom[ndom].wind_rho_max = zdom[ndom].rho_max = rho_max;
  zdom[ndom].zmax = zmax;
  zdom[ndom].zmin = zmin;

  zdom[ndom].rmax = rmax;
  zdom[ndom].rmin = rmin;
  zdom[ndom].wind_thetamin = zdom[ndom].wind_thetamax = 0.;

  /* Set up wind planes around the cells which in the wind.  This can be
   * smaller than the entire grid.*/

  zdom[ndom].windplane[0].x[0] = zdom[ndom].windplane[0].x[1] = 0;
  zdom[ndom].windplane[0].x[2] = zdom[ndom].zmin;

  zdom[ndom].windplane[0].lmn[0] = zdom[ndom].windplane[0].lmn[1] = 0;
  zdom[ndom].windplane[0].lmn[2] = 1;

  zdom[ndom].windplane[1].x[0] = zdom[ndom].windplane[0].x[1] = 0;
  zdom[ndom].windplane[1].x[2] = zdom[ndom].zmax;

  zdom[ndom].windplane[1].lmn[0] = zdom[ndom].windplane[0].lmn[1] = 0;
  zdom[ndom].windplane[1].lmn[2] = 1;


  return (0);
}


/* The next section calculates velocities.  We follow the hydro approach of
 * getting those velocities from the original grid.  This is really only
 * used for setting up the grid
 */


/* velocity for cylindrical coordinates only interpolates.  One has to
 * interpoalte for vgrad.
 *
 * This routine is dangerous because of the way it works, if one tries
 * to update wmain[].v. This is because it actually uses values in
 * wmain[v].  Ideally one would try to avoid maing such calls, but
 * ksl has not found a good way to do this.
 * */



/**********************************************************/
/**
 * @brief      The velocity at any position in an imported cylindrical
 * model
 *
 * @param [in] ndom   The domain of the imported model
 * @param [in] x   A position
 * @param [out] v   The velocity at x
 * @return     The speed at x
 *
 * @details
 * This routine interpolates on the values read in for the
 * imported model to give one a velocity
 *
 * ### Notes ###
 * In practice this routine is only used to initalize v in
 * wind structure.  This is consistent with the way velocities
 * are treated throughout Python.
 *
 **********************************************************/

double
velocity_cylindrical (ndom, x, v)
     int ndom;
     double *x, *v;
{
  int j;
  int nn;
  int nnn[4], nelem;
  double frac[4];
  double vv[3];
  double speed;
  coord_fraction (ndom, 0, x, nnn, frac, &nelem);
  for (j = 0; j < 3; j++)
  {
    vv[j] = 0;
    for (nn = 0; nn < nelem; nn++)
      vv[j] += wmain[zdom[ndom].nstart + nnn[nn]].v[j] * frac[nn];
  }

  speed = length (vv);

  /* Now copy the result into v, which is very necessary if refilling wmain.v */

  v[0] = vv[0];
  v[1] = vv[1];
  v[2] = vv[2];

  return (speed);
}




/**********************************************************/
/**
 * @brief      Get the density for an imported cylindrical model at x
 *
 * @param [in] ndom   The domain for the imported model
 * @param [in] x   A position
 * @return     The density in cgs units is returned
 *
 * @details
 * This routine finds rho from the imported model
 * at a position x.  The routine does not interpolate rho, but
 * simply locates the cell associated with x
 *
 * ### Notes ###
 * This routine is really only used to intialize rho in the
 * Plasma structure.  In reality, once the Plasma structure is
 * initialized we always interpolate within the plasma structure
 * and do not access the original data.
 *
 * This routine depends on the assumpition that x corresponds
 * to the mid-point of a cell.
 *
 **********************************************************/

double
rho_cylindrical (ndom, x)
     int ndom;
     double *x;
{
  double rho = 0;
  double r, z;
  int i, j, n;

  r = sqrt (x[0] * x[0] + x[1] * x[1]);
  z = fabs (x[2]);

  i = 0;
  while (z > imported_model.wind_z[i] && i < imported_model.mdim)
  {
    i++;
  }
  i--;

  j = 0;
  while (r > imported_model.wind_x[j] && j < imported_model.ndim)
  {
    j++;
  }
  j--;

  n = j * imported_model.mdim + i;

  rho = imported_model.mass_rho[n];


  return (rho);
}
