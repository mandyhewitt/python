/***********************************************************/
/** @file  setup_polarization.c
 * @author mh
 * @date   August, 2020
 *
 * @brief  Routines for reading in the settings for polarization
 *
 * Contains the function for reading in polarization
 * settings.
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "atomic.h"
#include "python.h"


/**********************************************************/
/**
 * @brief      Gets parameters from input file for
 * polarization settings.
 *
 * @return    0
 *
 * @details
 * Reads in polarization settings from file. These settings
 * are fully-documented in the Sphinx docs, and briefly here:
 *
 * ### Polarization.type ###
 * Sets whether or not to do polarization:
 *   'none': Off (default)
 *   'single': assumes incoming photons are unpolarized and
 *   singly electron scattered only. No resonant line scattering.
 *   'multiple': polarization status is initially unpolarised however
 *   Stokes q and Stokes u parameters are preserved on scattering
 *
 * ### Notes ###
 * TODO Add to Sphinx if required
 **********************************************************/
int
get_polarization_params (void) {

  char values[LINELENGTH], answer[LINELENGTH];

  rdpar_comment("Parameters for Polarization (if needed)");

  // ========== DEAL WITH BASIC POLARIZATION TYPE ==========
  strcpy(answer, "none");
  sprintf(values, "%d,%d,%d", POL_NONE, POL_SINGLE_SCATTER, POL_MULTIPLE_SCATTER);
  geo.polarization = rdchoice("Polarization.type(none,single,multiple)", values, answer);

  if (geo.polarization == POL_NONE)
    {
      Log("Polarization has been turned off \n");
    }

  if (geo.polarization == POL_SINGLE_SCATTER)
    {
      Log("Single Scattering has been turned on \n");
    }

  if (geo.polarization == POL_MULTIPLE_SCATTER)
  {
    Log("Multiple Scattering is not supported yet\n");
  }

  return (0);
}
