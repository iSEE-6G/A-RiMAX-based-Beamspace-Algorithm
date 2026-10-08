This is an extension of a RiMAX-based algorithm, namely PARAMING [R1]. PARAMING is a PARAmetric method for joint
angles and tiMING estimation that exploits the full spacetime structure of the system model. As a
result, PARAMING provides 3D sensing (AoA, AoD and ToA) information for each target and clutter component
with low complexity and high resolution, by leveraging model-based transformations. 

The extensions allows for a beamspace approach of PARAMING algorithm and is applicable for any generic radiation pattern whose array response is known. 


This script is part of the work done in iSEE-6G under University of Piraeus Research Center (UPRC). 
In this demo, the extension is set up to detect M = 3 targets. The process is as follows: 
   1) Generate beamspace CSI for M targets.
   2) Stack all beam pairs over subcarriers.
   3) Build a frequency-Hankel matrix.
   4) Estimate delays using the PARAMING/ESPRIT eigenvalue step.
   5) Remove delays by LS.
   6) Estimate 3D AoA/AoD by known-pattern matching.
   7) Estimate complex gains.
   8) Repeat over SNR and Monte Carlo trials.

Interested researchers are encouraged to extend the proposed framework. 


Corresponding Author: Harris K. Armeniakos, harmen@unipi.gr 



[R1] S. Naoumi, A. Bazzi, R. Bomfin and M. Chafii, "High-Resolution Sensing in Communication-Centric ISAC: Deep Learning and Parametric Methods," IEEE J. Sel. Areas Commun., vol. 44, pp. 2201-2216, 2026. 