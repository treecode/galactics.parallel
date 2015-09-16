/**********************************************************************/
/* Routine for minimizing ChiSquared using a downhill simplex method. */
/* (see Numerical Recipes for more information). Chi2 is determined   */
/* from fits of the rotation curve, velocity dispersion and surface   */
/* brightness data to results of a galaxy model produced via Kuijken  */
/* and Dubinski N-body realizations of the DFs.  The Chisqr also      */
/* incorporates a measure of the difference bet. target Bulge, Disk   */
/* and Halo masses with the calculated values.  Note that the target  */
/* Halo mass is calculated from  the observed rotation curve:         */
/* Mbulge + Mhalo + 1.213(Mdisk) = R(vcirc)^2/G = 2.74E11 at 30kpc    */
/**********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define DOSIMPLEX 1           /* If DOSIMPLEX=0, it does not minimize */
#define TIDALRAD 0.         /* Set to 0. if not to fit */
#define l 2                   /* Number of harmonics (2: testing, 10: real) */
#define DISKMASS 0.        /* Set to 0. if not to fit (typical=1.2e11) */
#define BULGEMASS 0.      /* Set to 0. if not to fit (typical=0.2e11) */
#define MAXITER 10
#define TOL 1.0e-4
#define RCUTOFF 30.
#define NPARM 10
#define BNPTCLS 80000        /* Number of bulge particles (80000) */
#define DNPTCLS 10000        /* Number of disk particles (10000) */
#define HNPTCLS 50000        /* Number of halo particles (50000) */
#define MSCALE 2.325e9  
/* (6.672e-11*1.989e30/1e5^2/3.086e19) */
/* -> v in 100km/s, d in kpc, M in units above (e.g. mdisk now *v^2, etc.) */

FILE *chi2file, *rcutfile;
int printout, fixdisk;

/*  function declarations */
void Galaxy(float psi0, float sig0, float q, float rcrk2, float ra, 
	    float mdisk, float rhob, float psicut, float sigb, float sfrac);
float ChiSqr(float psi0, float sig0, float q, float rcrk2, float ra, 
	     float mdisk, float rhob, float psicut, float sigb, float sfrac);
float calc(float *pars);

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[]) {
  
  float psi0, sig0, q, rcrk2, ra, mdisk, rhob, psicut, sigb, sfrac;
  float chi2[NPARM+1];
  float newchi2, newnewchi2, lastchi2min;
  float pars[NPARM+1][NPARM], fixpar[NPARM], h[NPARM], delta[NPARM];
  float sum[NPARM+1], midpt[NPARM];
  float newpars[NPARM], newnewpars[NPARM];
  int i, j, mini, maxi, done, increment, nfree;
  //////////////////////////////////////////////////////////////////////////

  // pars[i]     => initial guesses at parameters
  // fixpar[i]  => 0 to allow par to vary, 1 to fix
  // h[i]       => range over which to step (10% of value)
  // Note that mdisk is not used if disk mass is defined >0 above.
  pars[0][0] = 13.373;  fixpar[0] = 0;  h[0] = 1.3; // psi0   -30
  pars[0][1] =   7.075;  fixpar[1] = 0;  h[1] = 0.6; // sig0    3.35
  pars[0][2] =   3.544;  fixpar[2] = 0;  h[2] = 0.35; // q       1
  pars[0][3] =   1.;  fixpar[3] = 1;  h[3] = 0.04; // rcrk2   0.4
  pars[0][4] =   1.;  fixpar[4] = 1;  h[4] = 0.50; // ra      4.5
  pars[0][5] =  30.;  fixpar[5] = 0;  h[5] = 5.00; // mdisk   60 
  pars[0][6] =  9.947;  fixpar[6] = 0;  h[6] = 1.50; // rhob    6.4
  pars[0][7] = -19.125;  fixpar[7] = 0;  h[7] = 1.00; // psicut -18
  pars[0][8] =   2.02;  fixpar[8] = 0;  h[8] = 0.20; // sigb    2  
  pars[0][9] =   0.793;  fixpar[9] = 0;  h[9] = 0.08; // sfrac   0.75

    // MDISK:
  if (DISKMASS > 0.1) {
    pars[0][5] = 1.0 * DISKMASS / MSCALE;
    fixpar[5] = 1; h[5] = 0.00; 
  }    
  if (fixpar[5] == 1) fixdisk = 1;

    //////////////////////////////////////////////////////////////////////////

  done = 0;
  increment = 0;

  rcutfile = fopen("rcutoff.in", "w");
  fprintf(rcutfile, "%f\n", RCUTOFF);
  fclose(rcutfile);

  chi2file = fopen("chisqr.out", "w");
   
  printout = 1;

  chi2[0] = calc(pars[0]);

  if (DOSIMPLEX == 1) {
      
    for (i=1; i<(NPARM+1); i++) {
      for (j=0; j<NPARM; j++) pars[i][j] = pars[0][j];
      if (!fixpar[i-1]) {
	pars[i][i-1] = pars[i][i-1] + h[i-1];
	chi2[i] = calc(pars[i]);
      }
    }

    lastchi2min = 999;
    while (!done && (increment < MAXITER)) {

      increment++;
      mini = 0;
      maxi = 0;
    
      for (i=1; i<(NPARM+1); i++) {
	if (!fixpar[i-1]) {
	  if (chi2[i] > chi2[maxi]) maxi = i;
	  if (chi2[i] < chi2[mini]) mini = i;
	}
      }
      fprintf(chi2file, 
	      "\n\n%d:\n** Min/Max chi2: (%d)%.4f\t(%d)%.4f\n** Midpoint: ", 
	      increment, mini, chi2[mini], maxi, chi2[maxi]);
      fflush(chi2file);

      if (fabs(chi2[mini] - chi2[maxi]) < TOL) done = 1;
      
      for (j=0; j<NPARM; j++) {
	nfree = 1;
	sum[j] = pars[0][j];
	for (i=1; i<(NPARM+1); i++) {
	  if (!fixpar[i-1] && (i != maxi)) {
	    sum[j] += pars[i][j];
	    nfree++;
	  }
	}
	midpt[j] = sum[j]/nfree; 
	fprintf(chi2file, "%.3f\t", midpt[j]); fflush(chi2file);
      }


      fprintf(chi2file, "\nDelta: "); fflush(chi2file);
      for (i=0; i<NPARM; i++) {
	if (!fixpar[i]) {
	  delta[i] = pars[maxi][i] - midpt[i];
	}
	else  delta[i] = 0;
	newpars[i] = midpt[i] - delta[i];
	fprintf(chi2file, "%.0e\t", delta[i]);
      }
      fprintf(chi2file, "\n"); fflush(chi2file);
	
      newchi2 = calc(newpars);
	
      ////
      if (newchi2 < chi2[mini]) {
	for (i=0; i<NPARM; i++) {
	  newnewpars[i] = newpars[i] - delta[i];
	}
	newnewchi2 = calc(newnewpars);
	if (newnewchi2 < newchi2) {
	  for (j=0; j<NPARM; j++) {
	    pars[maxi][j] = newnewpars[j];
	  }
	  chi2[maxi] = newnewchi2;
	}
	else {
	  for (j=0; j<NPARM; j++) {
	    pars[maxi][j] = newpars[j];
	  }
	  chi2[maxi] = newchi2;
	}
	mini = maxi;
      }
      
      //
      else if (newchi2 > chi2[maxi]) {
	for (j=0; j<(NPARM+1); j++) {
	  if (j != mini) { 
	    for (i=0; i<NPARM; i++) {
	      if (!fixpar[i]) {
		pars[j][i] = 0.5*(pars[j][i]-pars[mini][i]) + pars[mini][i]; 
	      }
	    }
	    chi2[j] = calc(pars[j]);
	  }
	}
	// in this case, don't check tolerance limit since min doesn't change
	continue;
      }
	
      //
      else {
	for (i=0; i<NPARM; i++) {
	  pars[maxi][i] = newpars[i];
	}
	chi2[maxi] = newchi2;
	// in this case, don't check tolerance limit since min doesn't change
	continue;
      }
	
      // check if we are done 
      if (fabs(chi2[mini] - lastchi2min) < TOL) done = 1;
      lastchi2min = chi2[mini];

      fprintf(chi2file, "** Minimum chi2: %.4f\n", chi2[mini]);
      fflush(chi2file);

	
    }
    fprintf(chi2file, 
	    "\nFinal pars:\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n",
	    pars[mini][0], pars[mini][1], pars[mini][2], pars[mini][3], 
	    pars[mini][4], pars[mini][5], pars[mini][6], pars[mini][7]);
      
  }
  else {
    mini = 0;
  }

  fclose(chi2file);

}
//////////////////////////////////////////////////////////////////////////////
/////////////////////////////// FUNCTIONS ////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////


////////////////////////////// Make Galaxy ///////////////////////////////////

void Galaxy(float psi0, float sig0, float q, float rcrk2, float ra, 
	    float mdisk, float rhob, float psicut, float sigb, float sfrac) {

  FILE *infile;

  system("make clean");
  // write parameters to input file (see README for info):
  infile = fopen("in.dbh", "w");
  fprintf(infile, "y\n %f %f %f %f %f\n", psi0, sig0, q, rcrk2, ra);
  fprintf(infile, "y\n %f 5.4 40 0.3 1\n", mdisk);
  fprintf(infile, "y\n %f %f %f\n", rhob, psicut, sigb);
  fprintf(infile, ".5 2000\n");
  fprintf(infile, "%d\n", l);         /* use 10 */
  fclose(infile);

  //write in.bulge:
  infile = fopen("in.bulge", "w");
  fprintf(infile, "%f\t#streaming fraction\n", sfrac);
  fprintf(infile, "%d\t#number of particles\n", BNPTCLS);
  fprintf(infile, "-1\t#negative integer seed\n");
  fprintf(infile, "1\t#center the data 1=yes\n");
  fprintf(infile, "dbh.dat\t#harmonics file\n");
  fclose(infile);
  
  //write in.disk:
  infile = fopen("in.disk", "w");
  fprintf(infile, "%d\t#number of particles\n", DNPTCLS);
  fprintf(infile, "-1\t#random integer seed\n");
  fprintf(infile, "1\t#1=yes we want to center 0=no we don't\n");
  fprintf(infile, "dbh.dat\t#multipole expansion data file\n");
  fclose(infile);

  //write in.halo:
  infile = fopen("in.halo", "w");
  fprintf(infile, "0.5\t#streaming fraction\n");
  fprintf(infile, "%d\t#number of particles\n", HNPTCLS);
  fprintf(infile, "-1\t#negative integer seed for random number generator\n");
  fprintf(infile, "1\t#1=yes we want to center\n");
  fprintf(infile, "dbh.dat\t#multipole expansion data file\n");
  fclose(infile);

  // build the galaxy and output results to files:
  //system("rm -f vr.out mass_radius.out mass.out");
  system("rm -f fitvc.out fitsbp.out fitdisp.out");
  system("make galaxy");

  // /homes/titan/widrow/MEGA/DFs/GalactICS-exp/bin/
  system("bin/fitvc | tail -1 > fitvc.out");
  system("bin/vcirc > vr.out");
  system("bin/fitsbp > fitsbp.out");
  system("bin/fitdisp > fitdisp.out");
  //system("bin/mass_radius > mass_radius.out");
  system("bin/nfw > nfw.out");
    
  if (!DOSIMPLEX) {
    //system("bin/m_vs_r > mass.out");
    printf("Calculating axis ratios...\n");
    system("bin/axes > axes.out");
    printf("...writing to axes.out.\nDone!\n");
    //system("bin/los");
  }
  return;
}
//////////////////////////////////////////////////////////////////////////////

/////////////////////////// Minimize ChiSquared //////////////////////////////

float ChiSqr(float psi0, float sig0, float q, float rcrk2, float ra, 
	     float mdisk, float rhob, float psicut, float sigb, float sfrac) {

  FILE *fitvcfile, *fitsbpfile, *fitdispfile, *nfw, *tidalrfile;
  float chi2vc, massdisk, massbulge, masshalo, rhalo;
  float chi2sbp, mldisk, mlbulge, adisk, abulge;
  float cntr_min, rscntr_min;
  float chi2disp;
  float chi2fit, chi2;
  float sigmadisk, sigmabulge, sigmahalo, sigmatr, tidalr;
  float bfract, massrad, tmp, innermass, TOTALMASS, HALOMASS;
  int N;
    
  if (DISKMASS > 0) {
    // total mass is found from rotation curve at ~ 30kpc:
    TOTALMASS = RCUTOFF * 3.086e21 * 3.92e14 / (6.67e-08 * 1.989e33);
    HALOMASS = TOTALMASS - (1.213 * DISKMASS) - BULGEMASS;
  }
  else {
    HALOMASS = 0.;
    TOTALMASS = 0.;
  }

  fprintf(chi2file, 
	  "# psi0\tsig0\tq\trcrk2\tra\tmdisk\trhob\tpsicut\tsigb\tsfrac\n"); 
  fprintf(chi2file, 
	  "%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n",
	  psi0, sig0, q, rcrk2, ra, mdisk, rhob, psicut, sigb, sfrac);
  fflush(chi2file);

  Galaxy(psi0, sig0, q, rcrk2, ra, mdisk, rhob, psicut, sigb, sfrac);

  // read from input files
  fitvcfile = fopen("fitvc.out", "r");
  if (fitvcfile == NULL) {
    printf ("\nPROGRAM ABORTED:  fitvc.out does not exist.\n");
    exit(1);
  }
  fscanf(fitvcfile, "%f %f %f %f %f", 
	 &chi2vc, &massdisk, &massbulge, &masshalo, &rhalo);
  fclose(fitvcfile);
  //
  fitsbpfile = fopen("fitsbp.out", "r");
  if (fitsbpfile == NULL) {
    printf ("\nPROGRAM ABORTED:  fitsbp.out does not exist.\n");
    exit(1);
  }
  fscanf(fitsbpfile, "%f %f %f %f %f", 
	 &mldisk, &mlbulge, &adisk, &abulge, &chi2sbp);
  fclose(fitsbpfile);
  //
  fitdispfile = fopen("fitdisp.out", "r");
  if (fitdispfile == NULL) {
    printf ("\nPROGRAM ABORTED:  fitdisp.out does not exist.\n");
    exit(1);
  }
  fscanf(fitdispfile, "%f", &chi2disp);
  fclose(fitdispfile);
  //
  tidalrfile = fopen("tidalr.out", "r");
  if (tidalrfile == NULL) {
    printf ("\nPROGRAM ABORTED:  tidalr.out does not exist.\n");
    exit(1);
  }
  fscanf(tidalrfile, "%f", &tidalr);
  fclose(tidalrfile);
  if (TIDALRAD > 0) {
    sigmatr = fabs(tidalr-TIDALRAD)/(TIDALRAD) / 0.1;
  }
  else {
    sigmatr = 0.;
  }
    
  nfw = fopen("nfw.out", "r");
  if (nfw == NULL) {
    printf ("\nPROGRAM ABORTED:  nfw.out does not exist.\n");
    exit(1);
  }
  fscanf(nfw, "%f %f %f %f", &innermass, &cntr_min, &rscntr_min, &bfract);
  fclose(nfw);
  
  // determine "overall" chi2 by adding others in quadrature:
  chi2fit = sqrt(((chi2vc*chi2vc)+(chi2disp*chi2disp)+(chi2sbp*chi2sbp))/3);

  if (DISKMASS > 0.1)  
    sigmadisk = fabs(massdisk-DISKMASS)/DISKMASS / 0.02;
  else 
    sigmadisk = 0;

  if (BULGEMASS > 0.1) 
    sigmabulge = fabs(massbulge-BULGEMASS)/BULGEMASS / 0.02;
  else 
    sigmabulge = 0;

    // REMOVED DEPENDENCE ON HALOMASS:
    //sigmahalo = fabs(innermass-HALOMASS)/HALOMASS / 0.10;
  sigmahalo = 0;

  // How many "chisquareds" in minimization:
  N = 6;
  if (TIDALRAD == 0) N = N - 1;
  if ((DISKMASS == 0) || (fixdisk == 1)) N = N - 1;
  if (BULGEMASS == 0) N = N - 1;
       
  chi2 = sqrt(((chi2vc*chi2vc) + (chi2disp*chi2disp) + 
	       (chi2sbp*chi2sbp) + (sigmadisk*sigmadisk) + 
	       (sigmabulge*sigmabulge) + (sigmatr*sigmatr)) / (1.0 * N));

  // Write to galaxy file to watch evolution:
  if (printout) {
    //fprintf(chi2file, "%.4f\t", chi2); fflush(chi2file);
    fprintf(chi2file, 
	    "# Mdisk   Mbulge    Minner(%.0fkpc)   dM/L  bM/L\trhalo\tadisk\tabulge\n", RCUTOFF); 
    fprintf(chi2file, 
	    "%.3e %.3e %.3e       %.2f  %.2f\t%g\t%.2f\t%.2f\n", 
	    massdisk, massbulge, innermass, mldisk, mlbulge,
	    rhalo, adisk, abulge);
    fprintf(chi2file, 
	    "# DISKMASS=%.1e  BULGEMASS=%.1e  HALOMASS=%.1e  TOTALMASS=%.1e\n",
	    DISKMASS, BULGEMASS, HALOMASS, TOTALMASS);
    fflush(chi2file);  
    //fprintf(chi2file, "# RCUTOFF: %g kpc\n", RCUTOFF);
    if (TIDALRAD > 0) {
      fprintf(chi2file, 
	      "# TIDALRAD: %.2f kpc (%.0f requested)\n", tidalr, TIDALRAD);
      fflush(chi2file);  
    }
    else {
      fprintf(chi2file, 
	      "# TIDALRAD: %.2f kpc (not constrained)\n", tidalr);
      fflush(chi2file);  
    }
    //bfract = (massdisk + massbulge) / (massdisk + massbulge + masshalo);
    fprintf(chi2file, "# Baryon Fraction: %.2f\n", bfract);
      
    if (!DOSIMPLEX) {
      fprintf(chi2file, "# NFW: cntr_min = %.2f, rsmin*cntr_min = %.2f\n",
	      cntr_min, rscntr_min);
      fflush(chi2file);  
    }

    fprintf(chi2file, 
	    "# Chi2 ->\tvc: %.2f  \tsbp: %.2f\tbvel: %.2f\n \t\tmdisk: %.2f\tmbulge: %.2f\tmhalo(trunc): %.2f\n\t\trtidal: %.2f\n\t\tFIT: %.3f\tALL: %.3f\n", 
	    chi2vc, chi2sbp, chi2disp, 
	    sigmadisk, sigmabulge, sigmahalo, sigmatr, chi2fit, chi2);
    fflush(chi2file);  
  }

  return(chi2);
}
//////////////////////////////////////////////////////////////////////////////

float calc(float *pars) {
  return(ChiSqr(pars[0], pars[1], pars[2], pars[3], pars[4], pars[5], 
		pars[6], pars[7], pars[8], pars[9]));
}
//////////////////////////////////////////////////////////////////////////////
