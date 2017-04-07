#include "main.h"

float rad;
float Rc, zc;
extern void mysrand(const long seed);
int dotalk;

void gen_disk(int nptcl, int myseed, float *buffer, int verbose)
{
  dotalk = verbose;
	int i, j, k, nobj=10000;
	int seed= -123;
	int icofm=1;
	float q=0.9, rtrunc=5.0;
	float rhoran, massapp;
	float x, y, z, vx, vy, v, v2, R, vmax, vmax2;
	float phi, cph, sph, cth, sth, vR, vp, vz;
	float E, Lz, rad, rad2;
	float f0, frand, fmax, fmax1, vphimax1, psi, fnorm;
	float Fdisk();
	float ran1();
	float invu();
	float dr, rhomax1;
	float t, mass;
	float u1, v1, u1max, v1max;
	float xcm, ycm, zcm, vxcm, vycm, vzcm;
	float zip = 0.0, psi0;
	float rhomax=0.15, rhomin, rhotst;
	float rhoguess, zcon;
	float stream=0.5;
	float con, outdisk, drtrunc;
	float omega, kappa;
	float dv, v0;
	float broadcon=1.0;
	float vfacR, vfacz;
	float f1;
	float gr, gp, gz, g2;
	float vsigmax, vmean, vsig2;
	float gasdev();
	float diskdensf_(), sigr2_(), sigz2_();
	float FindMax(), FindMax1(), Fmax();
	float simpson(), massring();
	float velocityfactors_();
//	FILE *rhoout;
//	FILE *errfile;
	char harmfile[80];

//	errfile = fopen("errmsg.dat","w"); 
  if(verbose)
    fprintf(stderr,"The Disk\n");
	strcpy(harmfile,"dbh.dat");
  nobj = nptcl;
  seed = myseed;
  mysrand(myseed);
  icofm = 1;
  if (nobj < 0)
  {
    iquery("Enter the number of particles",&nobj);
    iquery("Enter negative integer seed",&seed);
    iquery("Center the simulation (0=no,1=yes)",&icofm);
    cquery("Enter harmonics file",harmfile);
  }
  

	r = (phase *) calloc(nobj,sizeof(phase));

	readdiskdf_(harmfile,&gparam);
	mdisk = gparam.mdisk; rdisk = gparam.rdisk; zdisk = gparam.zdisk;
	outdisk = gparam.outdisk; drtrunc = gparam.drtrunc;
	rd = 1.2*rdisk;
	zd = 1.2*zdisk;
	rtrunc = (outdisk + 2.0*drtrunc);
	diskmass = 4.0*M_PI*simpson(massring,0.0,rtrunc,128);
	mass   = diskmass/nobj;
	dr = rtrunc/100.;
	rhomax = 0.0;
//	rhoout = fopen("rhotst.dat","w");
	for(i=0; i<100; i++) {
	  R = i*dr;
	  z = 0.0;
	  rhoguess = exp(-R/rd);
	  rhotst = diskdensf_(&R,&z);
	  rhotst /= rhoguess;
	  if( rhotst > rhomax ) rhomax = rhotst;
//	  fprintf(rhoout,"%g %g\n",R, rhotst);
	}
	rhomin = 0.0;
	rhomax *=1.2;
#if 0
	{
		FILE *disksig;
		float rad, vsigR, vsigp;;

		disksig = fopen("disksig.dat","w");
		for(rad=0.1; rad<30.0; rad+=0.1) {
			omekap_(&rad,&omega,&kappa);
			vsigR = sqrt(sigr2_(&rad));
			vsigp = kappa/(2.0*omega)*vsigR;
			fprintf(disksig,"%g %g %g %g %g\n",rad,vsigp,vsigR,omega,kappa);
		}
		fclose(disksig);
	}
#endif

/*
	tstFdisk();
*/
  if (verbose)
    fprintf(stderr,"Calculating disk positions and velocities\n");
	vmean = 0;  vsig2 = 0;
	for(i=0, j=0, k=0; i<nobj;) {

	  u1 = -ran1(&seed);
	  v1 = 2.0*(ran1(&seed) - 0.5);
	  R = rd*invu(u1);
	  z = zd*atanh(v1);
	  
	  zcon = cosh(z/zd);
	  rhoguess = exp(-R/rd)/(zcon*zcon);
	  /* Guess at the approximate functional form of the density */
	  
	  rhotst = diskdensf_(&R,&z);
	  rhotst /= rhoguess;
	  
	  k++;
	  if( rhotst < rhomin )
	    continue;
	  
	  rhoran = (rhomax - rhomin)*ran1(&seed);
	  if( rhoran > rhotst )
	    continue;
	  phi = 2.0*M_PI*ran1(&seed);
	  x = R*cos(phi);
	  y = R*sin(phi);
	  omekap_(&R, &omega, &kappa);
	  vphimax = omega*R;	  
	  vsigR = sqrt(sigr2_(&R));
	  vsigp = kappa/(2.0*omega)*vsigR;
	  vsigz = sqrt(sigz2_(&R));
	  vsigmax = vsigR;
	  if( vsigp > vsigmax ) vsigmax = vsigp;
	  if( vsigz > vsigmax ) vsigmax = vsigz;


	  fmax = 1.1*FindMax(R,z,&vphimax);
/*
	  fmax1 = 1.1*FindMax1(R,z,&vphimax1);
	  fprintf(stderr,"fmax %g fmax1 %g vphimax %g vphimax1 %g\n",
			fmax, fmax1, vphimax, vphimax1);
*/
	  f0 = 0.0; frand = 1.0; /* dummy starters */
	  while( frand > f0 ) {
			  /*
		g2 = 999.;
	    while( g2 > 1.0) {
	      gr = 2.*(ran1(&seed) - 0.5);
	      gp = 2.*(ran1(&seed) - 0.5);
	      gz = 2.*(ran1(&seed) - 0.5);
	      g2 = (gr*gr + gp*gp + gz*gz);
	    }
		*/

		gr = 6.0*(ran1(&seed) - 0.5)*vsigR;
	    gp = 6.0*(ran1(&seed) - 0.5)*vsigp;
		gz = 6.0*(ran1(&seed) - 0.5)*vsigz;
/*
	    gp = (vphimax + 3.0*vsigp)*ran1(&seed) - vphimax;
*/
		/*
		gr = 3.0*vsigmax;
		gp = 3.0*vsigmax;
		gz = 3.0*vsigmax;
		*/
	    
	    vR = gr;
	    vp = vphimax + gp;
	    vz = gz;

	    f0 = Fdisk(vR, vp, vz, R, z);
	    frand = fmax*ran1(&seed);
#if 0
	    if( f0 > fmax ) {
	      float vpmax;
	      fprintf(errfile,"f0 > fmax at R=%g z=%g\nvr=%g vp=%g, vz=%g vphimax=%g f0=%g, fmax=%g\n", 
		      R,z, 
		      vR*broadcon/vsigR, 
		      (vp - vphimax)*broadcon/vsigp, 
		      vz*broadcon/vsigz, vphimax, f0, fmax);
	      fflush(errfile);
	    }
#endif
	    j++;
	  }
	  
	  velocityfactors_(&R, &z, &vfacR, &vfacz);
	  vphimax = vp - gp;
	  vphimax *= vfacR;
	  vp = vphimax + gp;
	  vz *= vfacz;
	  /*
	  fprintf(stdout,"%g %g\n",vp,frand);
	  */

	  cph = x/R; sph = y/R;
	  vx = vR*cph - vp*sph;
	  vy = vR*sph + vp*cph;
	  /*  vz stays the same */

    /* sanity check */
    if (isnan(x)) continue;
    if (isnan(y)) continue;
    if (isnan(z)) continue;
    if (isnan(vx)) continue;
    if (isnan(vy)) continue;
    if (isnan(vz)) continue;
    const float RMAX = 200;
    const float VMAX = 10;
    if (abs(x) > RMAX) continue;
    if (abs(y) > RMAX) continue;
    if (abs(z) > RMAX) continue;
    if (abs(vx) > VMAX) continue;
    if (abs(vy) > VMAX) continue;
    if (abs(vz) > VMAX) continue;

	  r[i].x = (float) x;
	  r[i].y = (float) y;
	  r[i].z = (float) z;
	  r[i].vx = (float)vx;
	  r[i].vy = (float)vy;
	  r[i].vz = (float)vz;
	  i++;
    if (verbose)
      if( i % 1000 == 0 ) {
        fprintf(stderr,".");
        fflush(stderr);
      }
  }
  if (verbose)
  {
    fprintf(stderr,"\n");
    fprintf(stderr,"number of density trials %d\n",k);
    fprintf(stderr,"number of velocity trials %d\n",j);
  }

  if( icofm ) {
    xcm = ycm =zcm = vxcm =vycm =vzcm = 0;
    for(i=0; i<nobj; i++) {
      xcm += r[i].x;
      ycm += r[i].y;
      zcm += r[i].z;
      vxcm += r[i].vx;
      vycm += r[i].vy;
      vzcm += r[i].vz;
    }
    xcm /= nobj; ycm /=nobj; zcm /= nobj;
    vxcm /= nobj; vycm /=nobj; vzcm /= nobj;

    for(i=0; i<nobj; i++) {
      r[i].x -= xcm;
      r[i].y -= ycm;
      r[i].z -= zcm;
      r[i].vx -= vxcm;
      r[i].vy -= vycm;
      r[i].vz -= vzcm;
    }
  }


  if (buffer == NULL)
  {
    t = 0.0;
#ifdef ASCII
    fprintf(stdout,"%d\n",nobj);
    for(i=0; i<nobj; i++) {
      fprintf(stdout,"% 15.7e % 15.7e % 15.7e % 15.7e % 15.7e % 15.7e % 15.7e\n", mass, r[i].x, r[i].y, r[i].z, r[i].vx, r[i].vy, r[i].vz);
    }
#else
    for(i=0; i<nobj; i++)  {
      fwrite(&mass,sizeof(float),1,stdout);
      fwrite(r+i,sizeof(phase),1,stdout);
    }
#endif
  }
  else
  {
    const int min_disk = 000000000;
    const int max_disk = 100000000;
    int NEL = 8;
    int i,pc;
    for (i = 0, pc= 0; i < nobj; i++, pc += NEL)
    {
      *((int*)&buffer[pc]) = min_disk + (i%(max_disk-min_disk));
      buffer[pc+1] = mass;
      buffer[pc+2] = r[i].x;
      buffer[pc+3] = r[i].y;
      buffer[pc+4] = r[i].z;
      buffer[pc+5] = r[i].vx;
      buffer[pc+6] = r[i].vy;
      buffer[pc+7] = r[i].vz;
    }
  }
}

float Fdisk(vr, vp, vz, R, z)
  float vr, vp, vz, R, z;
{
  float vr0, vp0, vz0, R0, z0;
  float diskdf5ez_();

  vr0 = vr; vp0 = vp; vz0 = vz; R0 = R; z0 = z;
  return diskdf5ez_(&vr0, &vp0, &vz0, &R0, &z0);
}

float invu(u)
  float u;
{
  /* invert the u function to find R */
  int i;
  double dr, rg, rnew, eps=1.0e-8;

  rg = 1.0; 
  dr = 1.0;
  eps = 1.0e-8;
  for(i=0; (i<20 && dr > eps); i++) {
    rnew = rg - (-(1.0 + rg)*exp(-rg) - u)
      /(rg*exp(-rg));
    dr = fabs(rnew - rg);
    rg = rnew;
  }
  if (dotalk)
    if( i == 20 ) fprintf(stderr,"Warning R did not converge\n");
  return (float) rg;
}

/* A approximate estimate of the local maximum of the distribution function */
float Fmax(R,z,vR,vz,vlimit, vpmax)
  float R, z, vR, vz;
  float *vpmax, vlimit;
{
  int i;
  float v, dv, fmax, vmax, v0;
  float Fdisk();

  dv = 2.0*vlimit/100;
  v0 = *vpmax - 50*dv;
  fmax = Fdisk(vR, v0, vz, R, z);
  for(i=0; i<100; i++) {
    float ftry, vtry;

    vtry = v0 + i*dv;
    ftry = Fdisk(vR, vtry, vz, R, z);
    if( ftry > fmax ) {
      fmax = ftry;
      vmax = vtry;
    }
  }
  *vpmax = vmax;

  return fmax;
}

float FindMax1(R,z,vpmax)

  float R, z, *vpmax;
{
  int upflag, flag;
  float dv, vpm, vpmold, v0, v1;
  float f0, f1, ftmp, fmid;
  float Fdisk();
  float zero;

  zero = 0.0;
  dv = 0.1*vsigp;
  vpm = *vpmax;
  v0 = vpm - dv;
  v1 = vpm + dv;

  f0 = Fdisk(0.0,v0,0.0,R,z);
  fmid = Fdisk(0.0,vpm,0.0,R,z);
  f1 = Fdisk(0.0,v1,0.0,R,z);
  if( fmid >= f0 && fmid >= f1 )
    return fmid;

  if( f0 > f1 ) {
    ftmp = f0;
    f0 = f1;
    f1 = ftmp;
    v1 = v0;
    dv = -dv;
  }
  vpm = v1;

  flag = 1;
  while( flag ) {

    dv *= 2.0;
    vpmold = vpm;
    vpm += dv;
    f0 = f1;
    f1 = Fdisk(0.0,vpm,0.0,R,z);
    flag = (f1 > f0);
  }
  *vpmax = vpmold;
  return f0;
}

float FindMax(R,z,vpmax)
  float R, z, *vpmax;
{
  float vpm;
  float v0, v1;
  float eps=0.01;
  float fmax, golden(), fgold();

  Rc = R; zc = z;

  vpm = *vpmax;
  v0 = *vpmax - 4.0*vsigp;
  v1 = *vpmax + 4.0*vsigp;
  fmax = -golden(v0,vpm,v1, fgold, eps, vpmax);

  return fmax;

}

float fgold(vp)
  float vp;
{
  float R, z;
  float zero;
  float Fdisk();

  zero = 0.0;
  R = Rc; z = zc;
  return( -Fdisk(0.0,vp,0.0,R,z) );
}

float massring(r)
  float r;
{
  float denz();
  float simpson();

  rad = r;
  return( rad*simpson(denz,0.0,5.0*zdisk,128) );
}

float denz(z)
  float z;
{
  float z0;
  float diskdensf_();

  z0 = z;
  return diskdensf_(&rad, &z0);
}

tstFdisk()
{
  float Fdisk(), FindMax();

  float R, z, vr, vp, vz, vphimax, fmax, omega, kappa;
  float v0, sigr2_(), vsigmax;
  FILE *fdiskfile;

  R=8.31; z = 0.2; vr = 0; vz = 0;

  omekap_(&R, &omega, &kappa);
  vphimax = omega*R;	  
  fmax = FindMax1(R,z,&vphimax);
  fdiskfile = fopen("fdisk.dat","w");
  fprintf(fdiskfile,"fmax %g vphimax %g vsigmax %g\n",fmax,vphimax,sqrt(sigr2_(&R)));
  vp = vphimax;
  vsigmax = sqrt(sigr2_(&R));
  for(v0=-3.0; v0<3.0; v0+=.002) {
    fprintf(fdiskfile,"%g %g %g %g %g\n",v0,v0+vphimax,
        Fdisk(0.0,v0+vphimax,0.0,R,0.0),
        Fdisk(0.0,v0+vphimax,0.0,R,0.01),
        Fdisk(0.0,v0+vphimax,0.0,R,0.02));
  }
  fclose(fdiskfile);
}
