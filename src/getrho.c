#include <stdio.h>
#include <stdlib.h>
#include <math.h>

main()
{
	float rhohalo();
	float radius, drho, rho0;

	readpotential();
	radius = 1.0;
	rho0 = rhohalo(radius,&drho);
	fprintf(stdout,"%g %g %g\n",radius,rho0, drho);
}
