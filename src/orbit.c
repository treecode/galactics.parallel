#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* This is a basic leapfrog orbit integration program that follows
   a test particle in the force field of a logarithmic potential -
   Phi = vc^2 log r where vc = 1
*/

float q1 = 0.9;
float q2 = 0.8;
float q1sq, q2sq;

main(argc,argv)
int argc;
char **argv;
{
	int i, n=20000;
	float x, y, z, vx, vy, vz, r, ax, ay, az, t, dt, ar, fr;
	float vxs, vys, vzs, binding_energy;
	float theta;
	float jx, jy, jz;
	float pot;
	float zlast;
	float m=1;
	float A, B, C, v0, q, psi0;
	FILE *phaseplot, *surfplot;
	void f(); 

/*
	phaseplot = fopen("phase.dat","w");
	surfplot = fopen("surfsec.dat","w");
*/
	readharmfile_("dbh.dat",&A,&B,&C,&v0,&q,&psi0);


/* Initial conditions */
	x = 5.;
	y = 7;
	z = 2.0;
	vx = 0.4;
	vy = -0.3;
	vz = -0.1;

	t = 0;
	dt = .025;
	
/* Calculate acceleration */

	r = sqrt(x*x + y*y); 
	force_(&r,&z,&ar,&az,&pot);
	ax = ar*x/r;
	ay = ar*y/r;

/* Initialize velocities by offsetting them by half a time step */
	vx += 0.5*ax*dt;
	vy += 0.5*ay*dt;
	vz += 0.5*az*dt;

	/*
	fprintf(stdout,"t x y z vx vy vz jx jy jz energy\n");
	*/

	zlast = z;
	for(i=0; i<n; i++) {

/* Update positions */
		x += vx*dt;
		y += vy*dt;
		z += vz*dt;
		r = sqrt(x*x + y*y);

/* Update accelerations */
		force_(&r,&z,&ar,&az,&pot);
		ax = ar*x/r;
		ay = ar*y/r;

/* Update velocities */
		vx += ax*dt;
		vy += ay*dt;
		vz += az*dt;
/* Synchronize velocities to determine total energy */
		vxs = vx - 0.5*ax*dt;
		vys = vy - 0.5*ay*dt;
		vzs = vz - 0.5*az*dt;

		binding_energy = 0.5*(vxs*vxs + vys*vys + vzs*vzs) + pot;

		jx = y*vzs - z*vys;
		jy = z*vxs - x*vzs;
		jz = x*vys - y*vxs;

		t += dt;
/*
		fwrite(&m,sizeof(float),1,stdout);
		fwrite(&x,sizeof(float),1,stdout);
		fwrite(&y,sizeof(float),1,stdout);
		fwrite(&z,sizeof(float),1,stdout);
		fwrite(&vxs,sizeof(float),1,stdout);
		fwrite(&vys,sizeof(float),1,stdout);
		fwrite(&vzs,sizeof(float),1,stdout);
*/
/*
		theta = atan2(y,x);
		if( theta < 0 ) theta += 2.0*M_PI;
		fprintf(stdout,"% 10.5e % 10.5e % 10.5e % 10.5e % 10.5e % 10.5e % 10.5e % 10.5e % 10.5e % 10.5e % 10.5e %10.5e %10.5e\n",t,x,y,z,vxs,vys,vzs,jx,jy,jz,binding_energy,theta,sqrt(vx*vx + vy*vy));
*/
		fprintf(stderr,"%g %g %g %g %g\n",t,x,y,z,binding_energy);
		zlast = z;
	}

}

void f(x,y,z,fx,fy,fz,pot)
float x, y, z, *fx, *fy, *fz, *pot;
{
	float r, m3, m2, m;
	float fr, epssq=.04;

	m2 = epssq + x*x + y*y/q1sq + z*z/q2sq;
	m = sqrt(m2);
	*fx = -x/m2;
	*fy = -y/m2/q1sq;
	*fz = -z/m2/q2sq;
	*pot = log(m);
}
