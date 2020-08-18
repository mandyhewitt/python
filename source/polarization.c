#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "atomic.h"
#include "python.h"


// define normalize function
// todo split the magnitude part out of this function

int normalize (double v1[], double v2[]) {
    double fMag;

    fMag = sqrt( pow(v1[0], 2) + pow(v1[1], 2) + pow(v1[2], 2));
    if (fMag == 0) {
        return 1;
    }
    else {
        v2[0] = v1[0] / fMag;
        v2[1] = v1[1] / fMag;
        v2[2] = v1[2] / fMag;
        return 0;
    }
}

double magnitude (double v1[]) {
    double Mag;
    Mag = sqrt( pow(v1[0], 2) + pow(v1[1], 2) + pow(v1[2], 2));
    return Mag;
}


int
polarization (p, pp)
    PhotPtr p;                              // the original wind photon structure
    PhotPtr pp;                             // the new extracted photon structure
{
    double *vect1;
    double *vect2;
    double adotb3;
    vect1 = p->orig;                         // vect1 now points to the address of the inbound photon
    vect2 = pp->lmn;                        // vect2 now points to the address of the outbound photon
    int i, j, mainloop, r1, c1, r2, c2, k, errno;
    double stokes_parameter_q_dimensionless, stokes_parameter_u_dimensionless, stokes_parameter_v_dimensionless;
    double thomson_scattering_angle, polarization_degree, intensity, diff;
    int loop;
    int ref_axis;                           // reference axis = 0 if z axis, 1 if x axis
    double prod1, prod2;
    double costheta, Iout;
    double magnitude_a,magnitude_b;
    double cosphi, sinphi;
    double phi_dot_result;
    double magnitude_r_obs;                  // magnitude of r axis in observer Frame of Reference
    double magnitude_r_out;                  // magnitude of r axis in scattering Frame of Reference
    double mag_mag;                          //
    double cross_mag;                        //
    double StokesVector[4][4] = {};          // Stokes vector
    double RotationMatrix[3][3] = {};        // Rotation matrix
    double StokesOut[4][1] = {};             // Stokes matrix scattered photon
    double StokesObs[3][1] = {};             // Stokes matrix scattered photon observer frame of reference
    double StokesIn[4][1] = {};              // Stokes Matrix for incoming photon
    double vectl[3];                         // l_scattering plane un-normalized vector
    double vectResult[3];                    // generic cross product result vector
    double vectrObsResult[3];                // unnormalized r_observer frame vector
    double vectlObsResult[3];                // unnormalized l_observer frame vector
    double vectNormal[3];                    // normalized r_scattering plane vector
    double vectlNormal[3];                   // normalized l_scattering plane vector
    double vectrObsNormal[3];                // normalized r_observer frame vector
    double vectlObsNormal[3];                // normalized l_observer frame vector*/
    double zaxis[3] = {0.0, 0.0, 1.0};       // z axis is the reference axis for our observer
    double xaxis[3] = {1.0, 0.0, 0.0};       // x axis is the alternative reference axis for our observer
    double vectxprod[3];                     // cross product r_obs and r_out to calculate sinphi
    double sinphi_numerator;                 // (r_obs x r_s).n_obs

    /* linearly polarized light is completely defined by the Stokes parameters Q and U (V=0) with the polarization
     * angle chi given by tan 2chi = U/Q and the angle gamma given by sin2gamma = V/I - but gamma is zero for
     * linearly polarized light */

    stokes_parameter_q_dimensionless = 0;    // Initially the photon is unpolarized
    stokes_parameter_u_dimensionless = 0;    // Initially the photon is unpolarized
    stokes_parameter_v_dimensionless = 0;    // Initially the photon is unpolarized
    thomson_scattering_angle = 0;            // Initialize Thomson scattering angle
    polarization_degree = 0;                 // Initially unpolarized degree
    adotb3 = 0.0;                            // dot product
    prod1 = 0;                               // Initialize dot product of x input position and x output position
    prod2 = 0;                               // Initialize dot product of y input position and y output position
    magnitude_a = 0;                         // Initialize magnitude of input vector
    magnitude_b = 0;                         // Initialize magnitude of output vector
    costheta = 0;                            // Initialize angle between the vectors
    cosphi = 0;                              // Initialize angle between scattered frame and observer frame of reference
    sinphi = 0;                              // Initialize
    intensity = p->w;                        // set intensity to photon weight
    loop = 100;                              // Loop counter for how many times to repeat the calculations
    mainloop = 0;                            // main program loop counter
    diff = 0;                                // difference between cos theta and 1
    errno = 1;                               // error detector
    Iout = 0;                                // Intensity of outgoing photon
    phi_dot_result = 0.0;

    /* Input position parameters will be in 3D cartesian coordinates on an x,y,z axis.
     * The first 3 coordinates entered represent the initial photon direction vector components.
     * The last 3 coordinates entered represent the final photon direction vector components. */

    adotb3 = (dot(vect1, vect2));

    magnitude_a = magnitude(vect1);                     // magnitude of input vector
    if (magnitude_a == 0)
    {
      Error ("Polarization: Input vector magnitude is zero \n");
      return (1);
    }

    magnitude_b = magnitude(vect2);                     //  magnitude of output vector
    if (magnitude_b == 0)
    {
        Error ("Polarization: Error in magnitude of output vector \n");
        return (1);
    }

    if (magnitude_a && magnitude_b)
    {
      //************************* LOCAL FRAME OF REFERENCE OF PHOTON **************************************
      // We have 2 non zero magnitudes so we can proceed to calculate our triad of axes and scattering angle
      // First calculate r scattering axis: normalized cross product between incoming vector and outgoing
      // vector
      //**************************** r axis local ********************************************************

      int m,n;
      ref_axis = 0;                                                   // initialize reference axis as 0 = z-axis
      n = memcmp( vect1, vect2, sizeof(*vect1) );

      if (n == 0)
      {
        Log ("Polarization: Input vector is the same as output vector.\n");    // input and output vectors are the same

        m = memcmp(vect1, zaxis, sizeof(*vect1));                   // check if they are the z axis

        if (m == 0) {
          Log("Polarization: Both vectors are z axis.\n");                  // input and output vectors are the z axis
          vectResult[0] = xaxis[0];                              // set outgoing scattering vector as x axis
          vectResult[1] = xaxis[1];                              // set outgoing scattering vector as x axis
          vectResult[2] = xaxis[2];                              // set outgoing scattering vector as x axis
          ref_axis = 1;                                          // reference axis flag is set to x
          Log("Polarization: Set reference axis as x axis.\n");              // input and output vectors are the x axis
        } else {
            vectResult[0] = zaxis[0];                               // set outgoing scattering vector as z axis
            vectResult[1] = zaxis[1];                               // set outgoing scattering vector as z axis
            vectResult[2] = zaxis[2];                               // set outgoing scattering vector as z axis
            Log ("Polarization: Set reference axis as z axis.\n");              // input and output vectors are the z axis
          }
        }
      else
        {
          cross(vect1, vect2, vectResult);
        }

      // now normalize the r axis in the scattering frame - new vector is called vectNormal
      if (normalize(vectResult,vectNormal) != 0) {
          // TODO: This is a very rare case and unlikely to happen in practice in the main code
          // but for completeness, this means that the input vector and the output vector are in the same
          // scattering plane.
          // - if r_s = r_obs i.e. the dotproduct of vect1 and vect2 = 0
          // then pick any vector perpendicular to vect2 e.g. parallel to or equal to r_obs.
          // Then print a warning to screen

          Error ("Polarization: Error in r scattering normalize function\n");
          Log ("Polarization ERROR Input and output vectors are anti-parallel\n");
          Log ("Polarization ERROR choose any perpendicular vector to output vector to be the r scattering axis\n");
          // TODO add this into the code or just skip this particular photon**************************!!!!!!!
          return (1);
        }

      //**************************** l axis local ********************************************************

      // Now we have r we can calculate l in the scattering plane
      // Calculate l scattering axis: normalized cross product between outgoing unit vector and r axis vector

      // todo is this the correct way round for this cross product??

      cross(vect2, vectResult, vectl);

      // now normalize the l axis in the scattering frame - new vector is called vectlNormal
      //todo the vector names aren't consistent - better to keep a theme going here

      if (normalize(vectl,vectlNormal) != 0) {
          Error ("Polarization: Error in l scattering normalize function\n");
          return (1);
        }

       //***********************OBSERVER FRAME OF REFERENCE***********************************************
       // Calculate r observer frame: normalized cross product between outgoing n vector and (0,0,1) ie z axis
       //**************************** r axis observer ********************************************************
       //todo add error checking

       if (ref_axis == 0) {
         cross(vect2, zaxis, vectrObsResult);
        } else {
          cross(vect2, xaxis, vectrObsResult);
        }

        // todo use python normalize command
        if (normalize(vectrObsResult,vectrObsNormal) != 0) {
            Log ("ERROR in Polarization r observer frame normalize function\n");
            Error ("Polarization:ERROR in r observer frame normalize function \n");
            return (1);                                   // TODO add error check on return
          }

        //**************************** l axis observer ********************************************************
        // Calculate l observer axis: normalized cross product between outgoing unit vector and r observer axis

        //todo add error checking

        cross(vect2, vectrObsNormal, vectlObsResult);

        if (normalize(vectlObsResult,vectlObsNormal) != 0) {
            Error ("Polarization: in l observer frame  normalize function\n");
            return (1);
          }

        //**************************SCATTERING ANGLE*****************************************************
        /* The cosine of the scattering angle between the input and output vectors
        * is equivalent to the dot product divided by the
        * multiplication of the input and output vector magnitudes: a · b = |a| × |b| × cos(θ) */

        // todo add error checking

        costheta = adotb3 / (magnitude_a * magnitude_b);

        // If the photon is passing straight through then angle = zero hence no point converting to degrees

        // TODO: check that this case is correct else fix the code so that costheta != 1.0

        diff = fabs(costheta) - 1.0000;

        if (fabs(diff) > 10e-8) {
          thomson_scattering_angle = (acos(costheta) * (180 / 3.14159265));  // convert angle to degrees
        }

        //**************************POLARIZATION DEGREE*****************************************************
        /* Now we have calculated the direction of propagation of the photon
        * and 2 vectors: inwards and outwards
        * Polarization in the plane can be calculated */

        // todo add error checking

        polarization_degree = (1 - (costheta * costheta)) / (1 + (costheta * costheta));

        //**************************STOKES VECTOR SETUP*****************************************************
        /* set up stokes matrix for unpolarized light*/
        /* Stokes Vector = [I, Q, U, V] nb we are not using v as it is for circularly polarized light only
        * Stokes Vector (Dimensionless) = Stokes Vector / Intensity
        * [1,q,u,v] = [1,0,0,0] for unpolarized light */
        //**************************INPUT STOKES VECTOR*****************************************************

        StokesIn[0][0] = intensity;

        // set up stokes vector

        StokesVector[0][0] = (0.75 * ((costheta * costheta) + 1));
        StokesVector[0][1] = (0.75 * ((costheta * costheta) - 1));
        StokesVector[1][0] = (0.75 * ((costheta * costheta) - 1));
        StokesVector[1][1] = (0.75 * ((costheta * costheta) + 1));
        StokesVector[2][2] = (1.5 * costheta);
        StokesVector[3][3] = (1.5 * costheta);

        // Set up some counters for the loops

        // element counters i j and k to go through the matrices

        i = 0;
        j = 0;
        k = 0;

        // row and column indicators for input and output matrices
        // todo calculate this rather than hard code it?

        r1 = 4;
        c1 = 4;
        r2 = 4;
        c2 = 1;

        //**************************OUTPUT STOKES VECTOR*****************************************************
        // Initializing all elements of StokesOut matrix to 0
        // todo is this even needed?

        for (i = 0; i < 4; ++i) {
            StokesOut[i][0] = 0;
        }

        // Multiplying matrices StokesArray and StokesIn and storing result in StokesOut matrix

        for (i = 0; i < r1; ++i) {
          for (j = 0; j < c2; ++j) {
            for (k = 0; k < c1; ++k) {
              StokesOut[i][j] += StokesVector[i][k] * StokesIn[k][j];
            }
          }
        }

        //**************************DIMENSIONLESS OUTPUT STOKES VECTOR*****************************************
        // Displaying the dimensionless stokes vector result

        Iout = StokesOut[0][0];
        for (i = 0; i < 4; ++i) {
          StokesOut[i][0] = (StokesOut[i][0] / Iout);
        }

        //**************************GET q,u,v *****************************************************************
        stokes_parameter_q_dimensionless = StokesOut[1][0];
        stokes_parameter_u_dimensionless = StokesOut[2][0];
        stokes_parameter_v_dimensionless = StokesOut[3][0];     // we don't use v here

        //**************************SANITY CHECK POLARIZATION DEGREE *****************************************
        // Calculate the polarization degree here - (a double check and should be equal to value calculated
        // above)

        polarization_degree = sqrt((stokes_parameter_q_dimensionless*stokes_parameter_q_dimensionless)
                + (stokes_parameter_u_dimensionless*stokes_parameter_u_dimensionless));

        //**************************PREPARE STOKES ROTATION MATRIX *****************************************

        //*****************CALCULATE PHI - THE CLOCKWISE ROTATION ANGLE BETWEEN SCATTER *********************
        // ************************PLANE AND MERIDIAN AND THE OBSERVER PLANE ********************************

        //cos phi is the dot product of r_obs and r_out divided by the magnitude of the 2 vectors

        //calculate magnitude_r_obs
        // todo logically this should come up above with the r obs calculations
        magnitude_r_obs = magnitude(vectrObsResult);

        if (magnitude_r_obs == 0)
        {
          Error ("Polarization: Null magnitude for r vector in observer frame of reference \n");
          return (1);
        }

        //calculate magnitude_r_out
        magnitude_r_out = magnitude(vectResult);

        if (magnitude_r_out == 0)
        {
          Error ("Polarization: Null magnitude for r vector out in photon frame of reference \n");
          return (1);
        }

        //calculate r_obs_dot_r_out
        // todo I think the normalised vectors would be fine here instead of dividing by the magnitudes

        phi_dot_result = dot(vectResult, vectrObsResult);
        mag_mag = magnitude_r_obs * magnitude_r_out;                    //the two magnitudes are already nonzero
        cosphi = (phi_dot_result / mag_mag);

        // calculate r_obs_cross_r_out
        cross(vectrObsResult, vectResult, vectxprod);

        // take scalar prod of vectxprod and direction to observer as both should be perp to dir of obs. par
        // or anti parallel

        cross_mag = magnitude(vectxprod);
        sinphi_numerator = dot(vectxprod,vect2);

        //todo vect 2 is not a unit vector so it is too big (but the sign is correct) calculate it!!
/*
            if (cross_mag == 0) {

                printf("ERROR null magnitude in cross product between vector r observer FoR and "
                               "vector r scattering frame of reference\n");
                return (1);                             // TODO what error codes should I be using?check python
            }

            printf("Magnitude of the cross product between vector r observer frame of reference and "
                           "vector r scattering frame of reference is: %f\n",cross_mag);
            */
        sinphi = (sinphi_numerator / (mag_mag*magnitude_b));

        // todo initialize rotation matrix

        // set the Rotation matrix

        RotationMatrix[0][0] = (1);
        RotationMatrix[1][1] = ((2 * (cosphi * cosphi)) - 1);
        RotationMatrix[1][2] = (2 * (sinphi * cosphi));
        RotationMatrix[2][1] = (-2 * sinphi * cosphi);
        RotationMatrix[2][2] = ((2 * cosphi * cosphi) - 1);

        // Displaying the Rotation matrix

        i = 0;
        j = 0;
        k = 0;

        r1 = 3;
        c1 = 3;
        r2 = 3;
        c2 = 1;
/*
            printf("\nRotation matrix:\n");

            for (i = 0; i < r1; ++i) {
                for (j = 0; j < c1; ++j) {
                    printf("%f  ", RotationMatrix[i][j]);
                    if (j == c1 - 1) {
                        printf("\n\n");
                    }
                }
            }
*/
        // Initializing all elements of StokesObs matrix to 0
        // todo is this even needed?

        for (i = 0; i < 4; ++i) {
          StokesObs[i][0] = 0;
        }

        // Displaying the result
/*
            printf("\nInitial Stokes Observer Matrix:\n");


            for (i = 0; i < 3; ++i) {
                    printf("%f  \n\n", StokesObs[i][0]);
            }
*/
        // todo - the matrix should be I Q U not i q u v as it is just now!!! I needs to be 1 I think.
        // (i)(1)
        // (q)(-1)
        // (u)(0)

        i = 0;
        j = 0;
        k = 0;

        r1 = 3;
        c1 = 3;
        r2 = 3;
        c2 = 1;

        for (i = 0; i < r1; ++i) {
          for (j = 0; j < c2; ++j) {
            for (k = 0; k < c1; ++k) {
              StokesObs[i][j] += RotationMatrix[i][k] * StokesOut[k][j];
            }
          }
        }
        // Displaying the result
/*
        printf("\nFinal Stokes Observer Matrix:\n");
        for (i = 0; i < 3; ++i) {
            printf("%f  \n\n", StokesObs[i][0]);
        }
*/
        pp->q = StokesObs[1][0];
        pp->u = StokesObs[2][0];
 //           printf("POLTEST: Observer q is %f  \n\n", pp->q);
 //           printf("POLTEST: Observer u is %f  \n\n", pp->u);
        if (pp->q < -1 || pp->q > 1) {
          Error("Polarization: Abnormal Stokes parameter q is %f \n", pp->q);
        }
        if (pp->u < -1 || pp->u > 1) {
          Error("Polarization: Abnormal Stokes parameter u is %f \n", pp->u);
        }
    } else {
      Error ("Polarization: A vector magnitude is zero. Unable to calculate Thomson scattering angle "
                           "or Polarization degree\n");
    }
    // end of main loop
//    fclose (fp);
//    if (verbosity == 10) {printf ("File created successfully\n");}
    return (0);
}