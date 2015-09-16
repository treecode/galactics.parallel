# galactics.parallel

A modified version of Galactics that can be used together with Bonsai to generate initial conditions using multiple processes 

2015-09-16  

This is a modified version of the Galactics package by Dubinski, Kuijken and Widrow
the copyright of this package belongs to them.

The changes we made allow us to include the galactics package as library into Bonsai 
and generate models in parallel on multiple nodes and then combine them into a single 
large model.


To compile Galactics:
 - cd src
 - make -f Makefile.[ifort/gnu]  use ifort if you have Intel Fortran compiler (3x faster than gfortran), otherwise use Makefile.gnu
 - cp libgengalaxy.a ../../  (to link it with Bonsai)



To generate Initial conditions:

	Take the following steps:
	- Setup your model as you would normally do in Galactics
	- Run make and wait untill the following files are generated:
		cordbh.dat dbh.dat freqdbh.dat mr.dat denspsibulge.dat denspsihalo.dat

	- Once the above files are generated you can interrupt Galactics.
	- Copy the above files to the folder from which you run Bonsai 
	- Modify/create the file 'component_numbers.txt'. This file is used to determine 
          the ratio's between the number of particles in the Halo, Bulge and Disk.

 	It should contain on the first line three numbers:
	nHalo nBulge nDisk

	Again these numbers will be converted to ratio's based on the actual number of particles that are requested using 
	the Bonsai parameter: --milkyway nParticlesPerNode .
	- Now you are ready to run Bonsai and generate the IC



When using this code make sure to cite the original authors:

K - Kuijken, K.; Dubinski, J. 1995  (1995MNRAS.277.1341K) 
Widrow & Dubinski 2005 		(2005ApJ...631..838W)  
Widrow, Pym and Dubinski 2008 ( 2008ApJ...679.1239W )


Evghenii Gaburov
Jeroen Bédorf

=============================================================


GalactICS Quickstart Manual

GalactICS is a package for generating realistic models of galaxies that
can optionally include a disk, bulge, dark halo and central supermassive 
blackhole.  The code works by running a suite of programmes that compute
a gravitational potential for a given mass model, compute a distribution
function for each component and then realize these components as
independent N-body files that can be joined together to make an N-body
model of a galaxy.  The computed gravitational potential can also be used
as a force field in test particle orbital calculations and other
applications.

The physical and mathematical details behind the construction of these
models is described in Widrow, L. and Dubinski, J. 2005, ApJ, in press


Building and Installation
-------------------------

This is super easy.

1. cd src
2. make
3. make install

This will build the code and copy the binaries to the directory ./bin.

To clean and rebuild just type (make clean; make; make install).

Building Models
--------------

This is best illustrated with the examples of M31 and the Milky Way in the
subdirectory ./models.  Here you will find a series of input parameter files 
for various programmes that describe the underlying models.  Check the 
README and README.parameters file in example of M31 for a detailed description 
of these parameters and how to build a galaxy.


--------------
John Dubinski
May-31-2005


