#include "main.h"

readpotential(dbhname)
char *dbhname;
{
	int i, j, k, nobj=10000;
	int seed= -123;
	float q=0.9, rtrunc=5.0;
	float rhoran, massapp;
	float x, y, z, vx, vy, v, v2, R, vmax, vmax2;
	float phi, cph, sph, cth, sth, vR, vp, vz;
	float E, Lz, rad, rad2;
	float f0, frand, fmax, psi;
	float dfhalo_(), halodens_(), pot_();
	float ran1();
	float dr, rhomax1;
	float t, mass;
	float u1, v1, u1max, v1max;
	float xcm, ycm, zcm, vxcm, vycm, vzcm;
	float zip = 0.0;
	float rhomax=0.15, rhomin, rhotst;
	float rhocur, rcyl;
	float stream=0.5, massfrac=1.0;
	float psic_halo_n, psi_n;
	float fcut, fac, ebh, psic_bh;
	float dfmax, emax, psibh;
	float b, c, d, int_const, sig2, aeff;
	float v02, vb2, coef1, coef2;
	int icofm=1;
	int jj;
	char harmfile[80];

	strcpy(harmfile,dbhname);
	/*	cquery("Enter harmonics file",harmfile); */

	readmassrad(); /* reads in mass and truncation radius of components */
	readharmfile_(harmfile,&gparam,&cparam,&bparam,&flags);
	readdenspsihalo_();

	c = gparam.c; v0 = gparam.v0; a = gparam.a; 
	cbulge = gparam.cbulge; v0bulge=gparam.v0bulge;
	psi0 = gparam.psi0; haloconst = gparam.haloconst;

	psic_halo = cparam.psic_halo;
	bhmass = bparam.bhmass;
	rtrunc = haloedge;
	mass = halomass/nobj;
	sig2 = -psi0;
	u1max = rtrunc;
	v1max = 0.5*M_PI;

	
}

float rhohalo(radius,drho)
float radius;
float *drho;
{
	float rcyl, z, rhocur, halodens_();
	float rho0, rho1, dr;
	/* drho is the derivative of rho wrt r */

	rcyl = radius;
	z = 0.0;

  	rhocur = halodens_(&rcyl,&z); 
	dr = 0.001*radius;
	rcyl = radius*1.0005;
  	rho1 = halodens_(&rcyl,&z); 
	rcyl = radius*0.9995;
  	rho0 = halodens_(&rcyl,&z); 
	*drho = (rho1 - rho0)/dr;

	return rhocur;
}
