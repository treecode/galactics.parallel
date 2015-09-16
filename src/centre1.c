/*
 * Routine to determine the geometrical centre of a spherical distribution
 * of stars.  The centre of mass of the entire system is given as the
 * initial guess.  The system is considered centred when the centre of mass
 * position within a sphere of given radius is less than .1% of that radius
 * If the given sphere has a radius of 9999 then the centre of mass of the
 * system is taken to be the system centre
 */

#include <stdio.h>
#include <math.h>

typedef struct {
	float x;
	float y;
	float z;
	float vx;
	float vy;
	float vz;
} phase;

extern float *rsq;
float *locmass;
float cmx, cmy, cmz;

centre(nobj,mass,x,y,z,rc,indx,nb)
int nobj, *indx, *nb;
float rc;
float *mass, *x, *y, *z;
{
	int i, *p, count, centred;
	float eps;
	int compin();

	locmass = mass;
	for(i=0; i<nobj; i++) 
		indx[i] = i;

	count = 0; centred = (rc == 9999.); /* if rc = 9999 then don't iterate */
	if( rc == 9999 ) {
		*nb = nobj;
		radii(nobj,x,y,z);
		for(i=0; i<nobj; i++) indx[i] = i;
		qsort(indx,nobj,sizeof(int),compin);
		return;
	}
	if( rc > 0 ) 
		spcofm(nobj,x,y,z,indx,9999.,nb);	/* find cofm of system */
	else
		rc = -rc;
	eps = .01*rc;
	while( !centred && count++ < 20) {
		radii(nobj,x,y,z);
		spcofm(nobj,x,y,z,indx,rc,nb);
		centred = ( fabs(cmx)<eps && fabs(cmy)<eps && fabs(cmz)<eps );
fprintf(stderr,"%g %g %g\n",cmx,cmy,cmz);
	}
	if( count == 11 || count == 10 )
		fprintf(stderr,"centre1:Warning didn't converge after 10 iterations\n");
skip:
	radii(nobj,x,y,z);
	*nb = nobj;
	for(i=0; i<nobj; i++) indx[i] = i;
	qsort(indx,nobj,sizeof(int),compin);
}

/* Determine the centre of mass of the system 
   and centre the objects on it */
  
spcofm(nobj,x,y,z,indx,rc,nb)
int nobj;
int *indx, *nb;
float rc;
float *x, *y, *z;
{
	int i, j, k, nbin();
	float mtot;
	
	*nb = nbin(nobj,rc,indx);
fprintf(stderr,"nb is %d\n",*nb);
	cmx = cmy = cmz = 0;
	mtot = 0;
	for(k=0; k< *nb; k++) {
		i = indx[k];
		cmx += locmass[i]*x[i];
		cmy += locmass[i]*y[i];
		cmz += locmass[i]*z[i];
		mtot += locmass[i];
	}
	cmx = cmx/mtot;
	cmy = cmy/mtot;
	cmz = cmz/mtot;

	for(j=0; j<nobj; j++) {
		x[j] -= cmx;
		y[j] -= cmy;
		z[j] -= cmz;
	}
}

radii(nobj,x,y,z)
int nobj;
float *x, *y, *z;
{
	register int k;
	for(k=0; k<nobj; k++) 
		rsq[k] = x[k]*x[k] + y[k]*y[k] + z[k]*z[k];
}

int nbin(nobj,rc,indx)
int nobj;
int *indx;
float rc;
{
	int i, j;
	float r2;

	j = 0;
	r2 = rc * rc;
	for(i=0; i<nobj; i++) {
		if( rsq[i] < r2 ) 
			indx[j++] = i;
	}
	return(j);
}

int compin(x,y)
int *x, *y;
{
	if(rsq[*x] > rsq[*y]) return(1);
	if(rsq[*x] < rsq[*y]) return(-1);
	return(0);
}
