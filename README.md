
# SSCSI-Coded-aperture-optimization 

This repository contains the optimization algorithm implemented in "Salazar et al., Coded Aperture Optimization in Spatial Spectral Compressive Spectral Imagers".
The original paper can be downloaded from https://ieeexplore.ieee.org/document/9034153. Please cite this paper if used. The cost function can be seen below ![Alt text](https://github.com/Edgar-Noita/SSCSI-Coded-aperture-optimization/blob/main/eq_22.png).


For details about the derivation process please refer to the original paper. A detailed explanation of the uploaded files can be seen below:

exe: Main file.

evaluation_2_blue_fin:  Function that executes the principal algorithm

ct_fn: Calculation of the first and second terms, EQ. (22)

ct_fn_ver: Calculation of the thrid term, EQ. (22)

ct_fn_2_opt_2: Calculation of the first and second terms, EQ. (22) for a single iteration

example: Implementation of the optimization for a given scenario. Plots the evolution of the cos function and also the patterns after each iteration

my_bool: Creates an array of matrices with boolean complementary codes (see restiction on Eq. (22)).

The figure below shows the evolution of the patterns over the iterations:

![Alt text](https://github.com/Edgar-Noita/SSCSI-Coded-aperture-optimization/blob/main/pic1.png)



