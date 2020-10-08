#!/bin/bash

#---PRINTS OUT THE HARTEE ENERGIES OF EACH STATE CORRESPONDING TO A GIVEN ALPHA VALUE---

for i in `seq 1 10`

do

	grep -i "MCSCF STATE" gs-irc-$i.out > copy1-$i.dat

	sed --in-place /Dipole/d ./copy1-$i.dat

	cut -c 35-58 copy1-$i.dat > energy-$i.dat

	rm copy1-$i.dat

done

#----------------------------------------------------------------------------------

#---PRINTS THE VARIATION IN THE ENERGY OF EACH STATE AS A FUNCTION OF ALPHA INTO DIFFERENT FILES--

for i in `seq 1 15`

do

	for j in `seq 1 10`

	do

		sed -n $[$i]p energy-$j.dat > hartree

		paste hartree >> gs-hartree-$i.dat

	done

done

#----------------------------------------------------------------------------------
