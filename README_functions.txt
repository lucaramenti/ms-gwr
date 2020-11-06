In this file we will briefly explain how to use the functions used in file functions.R

#bw_cv
This function is used to compute the optimal bandwidth, using LOO cross-validation or k-fold cross-validation, considering both bandwidths to be equal (this function can be easily generalized in order to set different bandwidths for site and event coordinates). At first an interval of possible bandwidths has to be selected, setting the minimum and maximum bandwidth (bw_min, bw_max) and a step (step), then the number of folds (f) has to be selected. If LOO cross-validation has to be carried out, parameter f has to be set qual to the number of data in the dataset. The next step is to decide if we want to use either SEC_only_calibration or ESC_only_calibration (func), explicitating the desired method in parameter "method" selecting either "sec" or "esc". The following required parameters (Xc, Xe, Xs, y, intercept, utm_ev_sp, utm_st_sp) are the input variables used for the execution of SEC_only_calibration or ESC_only_calibration. This function returns a list containing the optimal bandwidth and the corresponding cross-validated sum of squares.

#bw_gcv
This function is used to compute the optimal bandwidth, over a set of possible bandwidths, using the generalized cross-validation criterion (GCV) in its classical version. The input parameters are the same as in bw_cv.

#bw_gcv_mei
This function is used to compute the optimal bandwidth, over a set of possible bandwidths, using the generalized cross-validation criterion (GCV) in the version proposed by Mei. The input parameters are the same as in bw_gcv.

#gcv_mei_only_one
This function is used to compute the optimal bandwidth using the generalized cross-validation criterion (GCV) in the version proposed by Mei. Unlike the previous function, in this function we set the bandwidths (which can differ from each other) as an input parameter: gcv_mei_only_one can be cycled over different values of bwe and bws.

#SEC_general
This function is used to carry out a complete calibration and compute all regression coefficients over a given grid using the estimation order SEC. Xc, Xe, Xs, y are matrices (or vectors) containing all calibration covariates: the stationary ones, the spatially varying ones and the response variable. Then we have to decide whether the intercept (intercept) has to be constant ("c"), event-dependent ("e") or site-dependent ("s"). The next step is to set the desired bandwidths (bwe, bws), selected using one of the previous functions. Finally we have to give as an input the site- and event-coordinates (utm_st_sp, utm_ev_sp) and the grid over which the whole calibration has to be carried out. This function returns a list containing all the estimated regression coefficients and the matrices He, Hs and B, which can be useful to compute the hat matrix.

#ESC_general
This function is used to carry out a complete calibration and compute all regression coefficients over a given grid using the estimation order ESC. Input and output are the same as in SEC_general.

#SCE_general
This function is used to carry out a complete calibration and compute all regression coefficients over a given grid using the estimation order SCE. Input and output are the same as in SEC_general.

#ECS_general
This function is used to carry out a complete calibration and compute all regression coefficients over a given grid using the estimation order ECS. Input and output are the same as in SEC_general.

#CSE_general
This function is used to carry out a complete calibration and compute all regression coefficients over a given grid using the estimation order CSE. Input and output are the same as in SEC_general.

#CES_general
This function is used to carry out a complete calibration and compute all regression coefficients over a given grid using the estimation order CES. Input and output are the same as in SEC_general.

#emgwr_prediction_points
This function is used to make predictions over a whole set of points, using saved result from previous calibrations. Xc, Xe, Xs, y are data used for calibration and emgwr is output of a calibration method, which is needed in order to retrieve the resulting hat matrix. In this case parameters delta1 and delta2 have been saved in advance after carrying out the calibration, but it is also possible to compute them inside the function (see emgwr_prediction_points_no_param). bwe, bws, utm_ev_sp, utm_st_sp are the other parameters used to calibrate emgwr. The next step is to set the prediction covariates (pc, pe, ps) with their respective coordinates (pcoords, which is passed as a n x 4 matrix). Finally we have to set a confidence level alfa. This function returns the predicted values, together with the associated variance, the estimated regression coefficients and other parameters.

#emgwr_prediction_points_no_param
This function is used to make predictions over a whole set of points, using saved result from previous calibrations, with the only difference that in this case delta1 and delta2 are computed within the function and are not passed as an input.

#ESC_only_calibration
This function is used to carry out a complete calibration without computing all regression coefficients, using the estimation order ESC. The input parameters are the same as in SEC_general. This function returns the matrices which are needed to compute the resulting hat matrix.

#ESC_only_constant_intercept_calibration
This function is used to carry out a complete calibration without computing all regression coefficients, using the estimation order ESC, under the particular circumstance in which only the intercept is spatially staionary. The input parameters are the same as in ESC_only_calibration, except for Xc, which is not present. This function returns the matrices which are needed to compute the resulting hat matrix.

#ESC_no_intercept_calibration
This function is used to carry out a complete calibration without computing all regression coefficients, using the estimation order ESC, under the particular circumstance in which there is no intercept. The input parameters are the same as in ESC_only_calibration. This function returns the matrices which are needed to compute the resulting hat matrix.

#ESC_grid_creation
This function is used to compute the regression coefficients over a grid starting from saved results of a previous calibration, using the estimation order ESC. Xc, Xe, Xs, y, intercept, bwe, bws, utm_ev_sp, utm_st_sp are the input parameters used to calibrate the object emgwr, while grid is the grid over which we want to compute the regression coefficients.

#SEC_only_calibration
This function is used to carry out a complete calibration without computing all regression coefficients, using the estimation order SEC. The input parameters are the same as in SEC_general. This function returns the matrices which are needed to compute the resulting hat matrix.

#SEC_only_constant_intercept_calibration
This function is used to carry out a complete calibration without computing all regression coefficients, using the estimation order SEC, under the particular circumstance in which only the intercept is spatially staionary. The input parameters are the same as in SEC_only_calibration, except for Xc, which is not present. This function returns the matrices which are needed to compute the resulting hat matrix.

#SEC_no_intercept_calibration
This function is used to carry out a complete calibration without computing all regression coefficients, using the estimation order SEC, under the particular circumstance in which there is no intercept. The input parameters are the same as in SEC_only_calibration. This function returns the matrices which are needed to compute the resulting hat matrix.

#SEC_grid_creation
This function is used to compute the regression coefficients over a grid starting from saved results of a previous calibration, using the estimation order SEC. Xc, Xe, Xs, y, intercept, bwe, bws, utm_ev_sp, utm_st_sp are the input parameters used to calibrate the object emgwr, while grid is the grid over which we want to compute the regression coefficients.

#gauss_kernel
This function returns a vector of weights for a Gaussian kernel, setting as input a vector of distances (d) and a bandwidth (h).

#mixed_SC_calibration
This function is used to carry out the classical MGWR, where the spatial non-stationarity depends only on one set of spatial coordinates.

#mixed_SC_no_intercept_calibration
This function is used to carry out the classical MGWR, where the spatial non-stationarity depends only on one set of spatial coordinates, in the particular case where there is no intercept.