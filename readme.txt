SwARM: Software for Absolute and Relative Motion (Satellite Theories)

Author: Bharat Mahajan
https://github.com/princemahajan

SwARM is a second-order analytic propagator for absolute and relative motion of satellites in the presence of gravitational harmonic perturbations. This code is part of the supplementary material to the article titled "State Transition Matrix for Satellite Relative Motion in the Presence of Gravitational Perturbations" authored by Bharat Mahajan, Srinivas R. Vadali, and Kyle T. Alfriend in the AIAA Journal of Guidance, Control, and Dynamics, submitted in December 2018.

Copyright 2018 Bharat Mahajan

It is a Jagerware. That is, it is free software provided "as is" without warranty of any kind, express or implied. In no event shall the authors or copyright holders be liable for  anything related to this software. The software is free to use, modify and/or redistribute as far as the activities are non-commercial. In any case, you should cite the original author and the publications (given below). If you find it useful, then you owe the copyright holder a Jagerbomb!

--------------------------
Installation and Execution
--------------------------

Install the SwARM.mlapp app in MATLAB. It is created and tested on MATLAB R2017b version. Provide all the required filenames with complete paths. Use '\\' as the path separator. 

For GMAT, use correct GMAT working directory path in "GMAT Working Dir" textbox. SwARM needs it to be able to find the GMAT propagation result files (ChiefStates.txt, DepStates.txt, ChiefOE.txt, DepOE.txt), which are created by the GMAT script, in the GMAT_Working_Dir\\output folder. 

"MAT Results File" is the name of a file the user must provide, in which the propagation results are stored. The values of the "Gravitational Parameter" and "Radius" must match the gravity model provided by the "Gravity Coefficient File."

If the code gives an error due to MEX files, then uncheck the "Use J2 MEX Code" checkbox. In case of some other error, try restarting the application. As a last resort, unpack the source code from the SwARM.mlapp app file and run the codes manually by execuing the file 'SwARMmain.m' in MATLAB. For more details, see the help menu in the app. 

--------------------------
Known Issues
--------------------------

1. If the tesseral harmonics are enabled and the reference orbit (Chief's orbit) for formation flying has zero eccentricity or the inclination, then relative motion analytic propagation has singularity issues. However, there are no issues with any zonal harmonic.

----------------------------------------------------
All previous publications on which the code is based
----------------------------------------------------

[1] Mahajan, B., "Absolute and Relative Motion Satellite Theories for Zonal and Tesseral Gravitational Harmonics," PhD Dissertation, Texas A&M University, 2018.

[2] Mahajan, B., Vadali, S.R. & Alfriend, K.T. Celest Mech Dyn Astr (2018) 130: 25. https://doi.org/10.1007/s10569-018-9818-8

[3] Bharat Mahajan, S. R. Vadali, and K. T. Alfriend, “Analytic solution of perturbed relative motion with zonal and tesseral harmonics,” In J. W. McMahon, Y. Guo, F. A. Leve and J. A. Sims (Eds.), Spaceflight Mechanics 2017: Proceedings of the 27th AAS/AIAA Space Flight Mechanics Meeting, held February 5-9, 2017 San Antonio, Texas, U.S.A. (pp. 1117-1134). [AAS 17-475] (Advances in the Astronautical Sciences; Vol. 160). San Diego, California: Published for the American Astronautical Society by Univelt, Inc. (2017).


[4] Bharat Mahajan, S. R. Vadali, and K. T. Alfriend, “Exact normalization of the tesseral harmonics,” In J. W. McMahon, Y. Guo, F. A. Leve and J. A. Sims (Eds.), Spaceflight Mechanics 2017: Proceedings of the 27th AAS/AIAA Space Flight Mechanics Meeting, held February 5-9, 2017 San Antonio, Texas, U.S.A. (pp. 2569-2588). [AAS 17-473] (Advances in the Astronautical Sciences; Vol. 160). San Diego, California: Published for the American Astronautical Society by Univelt, Inc. (2017)


[5] Bharat Mahajan, S. R. Vadali, and K. T. Alfriend, “Analytic solution for satellite relative motion: The complete zonal gravitational problem,” In R. Zanetti, R. P. Russell, M. T. Ozimek and A. L. Bowes (Eds.), Spaceflight Mechanics 2016: Proceedings of the 26th AAS/AIAA Space Flight Mechanics Meeting, held February 14-18, 2016, Napa, California, U.S.A. (pp. 3325-3348). [AAS 16-262] (Advances in the Astronautical Sciences; Vol. 158). San Diego, California: Published for the American Astronautical Society by Univelt, Inc. (2016)


[6] Bharat Mahajan, S. R. Vadali, and K. T. Alfriend, “Analytic solution for satellite relative motion with zonal gravity perturbations,” In M. Majji, J. D. Turner, G. G. Wawrzyniak, and W. T. Cerven (Eds.), Astrodynamics 2015: Proceedings of the AAS/AIAA Astrodynamics Specialist Conference held August 9–13, 2015, Vail, Colorado, U.S.A. (pp. 3583-3598). [AAS 15-705] (Advances in Astronautical Sciences; Vol. 156). San Diego, California: Published for the American Astronautical Society by Univelt, Inc. (2016).


