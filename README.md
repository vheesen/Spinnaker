#Spinnaker

============

SPectral INdex Numeral Analysis of K(c)osmic-ray Electron Radio-emission

SPINNAKER is a 1D cosmic-ray transport code that numerically solves equations for pure advection and diffusion and calculates synthetic non-thermal radio continuum profiles.

Homepage: http://www.heesens.de/Spinnaker.html

Reference: Heesen, Dettmar, Krause, Beck and Stein 2016, MNRAS, 458, 322

Papers that use the code:
- Heesen et al. 2018a, MNRAS, 474, 5049. Application to a radio galaxy. Extension to a quasi-1D model with varying radius and velocity, including adiabatic losses.
- Heesen et al. 2018b, MNRAS, 476, 1580. Application to 12 late-type edge-on galaxies with constant wind speed.
- Heesen et al. 2018c, MNRAS, 476, 1756. Application to IC 10, a nearby starburst dwarf irregular galaxy, with accelerating wind.
- Schmidt et al. 2018, in prep. Application to CHANG-ES galaxies NGC 891 and 4565.
- Heald et al. 2018, in prep. Application to LOFAR and CHANG-ES data of the galaxy NGC 5775.
- Stein et al. 2018, in prep. Application to CHANG-ES data of the galaxies NGC 4013, 4217, and 4666.

PhD theses that use the code:
- Stein 2017, Ruhr-University Bochum, https://hss-opus.ub.ruhr-uni-bochum.de/opus4/frontdoor/index/index/start/3/rows/10/sortfield/score/sortorder/desc/searchtype/simple/query/Stein/docId/5556

============

*Documentation*

In order to run SPINNAKER you will have to compile it by running compile.sh, which created the file spectral.x. This is an executalbe that can be run by ./spectral.x. You may have to do 'chmod u+x compile.sh' and 'chmod u+x spectral.x'. Some preliminary documentation can be found in the 'doc' folder.

There is now also an interactive Version of Spinnaker available programmed in Python 'Spinteractive'. This was contributed by Arpad Miskolczi. All you need is spinteractive3.py and parameters-template. You may need to install the following packages: python-tk,  numpy, matplotlib. The rest should be installed by default. It is written for Python 2.7, not 3.

The input files that you wish to fit your parameters to needs to have the following format:

z; I; Delta_I
....

where z refers to the distance from the center in parsec, I is the intensity at that point and Delta_I the uncertainty of I. The unit of I and Delta_I does not matter because it is normalized to 1 anyway. At least two input files are needed at two different frequencies.

Futher documentation about Spinteractive is forthcoming.
