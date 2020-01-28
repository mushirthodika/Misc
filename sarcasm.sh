#!/bin/bash

echo ""
echo "*************************************************************"
echo ""
echo " SARCasM:   Simulating AbsoRption spectra of Complex Molecules "
echo ""
echo "    Written for my one and only true love: Barbara Marchetti   "
echo ""
echo "*************************************************************"
echo ""
SPECDIR=/users/chxtk/SArCASM/bin

rm -rf struct-*.xyz

cp -rf $SPECDIR/initqp_input .

echo -n "Name of file containing input geometry?    "

        read GEOM

                $SPECDIR/xyz2nx < $GEOM

echo ""

echo -n "Number of Atoms?    "

        read ATOMS

                TOT=$(($ATOMS + 1))

                sed -i "s/numat = 3/numat = $ATOMS/g" initqp_input

echo ""

echo -n "Number of Wigner Points?     "

        read POINTS

                sed -i "s/npoints = 10/npoints = $POINTS/g" initqp_input

echo ""

echo -n "Number of Gaussian log file containing Frequencies?    "

        read LOGFILE

                sed -i "s/file_nmodes = h2o.log/file_nmodes = $LOGFILE/g" initqp_input
echo ""

echo -n "Level of theory?    "

        read LEVEL

echo ""

echo -n "Number of absorbing/emitting states?   "

        read STATES


$SPECDIR/initcond.pl > initcond.log

rm -rf initcond.log
rm -rf DEBUG
rm -rf TEMP
rm -rf initgp
rm -rf initqp_input

mv final_output geometries

for i in $(seq 1 $POINTS)

        do

                grep -w -A$TOT "Initial condition =     $i" geometries | tail -n$ATOMS > struct-$i.xyz

                                grep -w -A$TOT "Initial condition =    $i" geometries | tail -n$ATOMS >> struct-$i.xyz

                                grep -w -A$TOT "Initial condition =   $i" geometries | tail -n$ATOMS >> struct-$i.xyz

                                awk '{print $1"\t"$3"\t"$4"\t"$5}' struct-$i.xyz > structure-$i.xyz

                                rm -rf struct-$i.xyz

                        done

                for i in $(seq 1 $POINTS)

                        do

                        echo "# td(nstate=$STATES) $LEVEL units(bohr)" > top

                        echo " " >> top

                        echo "title " >> top

                        echo " " >> top
echo "0 1" >> top

                        cat top structure-$i.xyz  > point-$i.com

                        echo " " >> point-$i.com

                        done

                rm -rf structure-*.xyz

                rm -rf geometries

                echo ""
                echo "😗 BACIO BACIO PRINCIPESSA 😗"
