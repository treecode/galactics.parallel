#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define GCON 6.349363592e-2 /* (2.0*pi)^(-3/2) */

/* Generate a distribution which follows an oblate logarithmic potential */

typedef struct {
  float x, y, z, vx, vy, vz;
}  phase;

typedef struct {
  float c;
  float v0;
  float a;
  float cbulge;
  float v0bulge;
  float abulge;
  float psi0;
  float haloconst;
  float bulgeconst;
  float mdisk;
  float rdisk;
  float zdisk;
  float outdisk;
  float drtrunc;
  float potcor;
} galaxy;

typedef struct {
  float psic_bulge;
  float psic_halo;
} cutoff;

cutoff cparam;

typedef struct {
  float bhmass;
  float aeff;
} blackhole;

blackhole bparam;

typedef struct {
  int idiskflag;
  int ibulgeflag;
  int ihaloflag;
  int ibhflag;
} components;

components flags;

phase *r;
galaxy gparam;


/* Disk Constants */
float mdisk, rdisk, zdisk, outdisk, drtrunc;
float rd, zd;
float vsigR, vsigp, vsigz, vphimax;
float diskmass, diskedge;

/* Bulge Constants */
float cbulge, v0bulge, abulge, bulgeconst;
float bulgemass, bulgeedge;
float psic_bulge;

/* Halo Constants */
float c, v0, a, haloconst, v02, psi0,psihalocut;
float halomass, haloedge;
float psic_halo;

/* Blackhole Constants */
float bhmass, bhsoft;
