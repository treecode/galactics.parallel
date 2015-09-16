#Build galaxy models
BIN=../../bin

galaxy: disk bulge halo
	cat disk bulge halo > galaxy

galaxy_ascii: galaxy
	$(BIN)/toascii galaxy > galaxy_ascii.dat

potential: dbh.dat

disk: freqdbh.dat cordbh.dat dbh.dat in.disk
	$(BIN)/gendisk < in.disk > disk

bulge: dbh.dat in.bulge
	$(BIN)/genbulge < in.bulge > bulge 

halo: dbh.dat in.halo
	$(BIN)/genhalo < in.halo > halo

dbh.dat: in.dbh
	$(BIN)/dbh < in.dbh

freqdbh.dat: dbh.dat h.dat b.dat
	$(BIN)/getfreqs

cordbh.dat: dbh.dat freqdbh.dat in.diskdf
	$(BIN)/diskdf < in.diskdf

vr.dat: dbh.dat h.dat b.dat
	$(BIN)/vcirc 1 1 > vr.dat

clean:
	rm -f disk bulge halo galaxy

veryclean:
	rm -f *.dat disk bulge halo galaxy tidalr.out
