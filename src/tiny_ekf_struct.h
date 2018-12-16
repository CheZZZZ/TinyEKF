/*
 * tiny_ekf_struct.h: common data structure for TinyEKF
 *
 * You should #include this file after using #define for N (states) and M
*  (observations)
 *
 * Copyright (C) 2016 Simon D. Levy
 *
 * MIT License
 */

typedef struct {

    int nn;          /* number of state values */
    int ne;          /* number of state values */
    int m;          /* number of observables */

    double x[NNsta];    /*nominal state vector */
    double dx[NEsta];   /*error-state vector*/
    double qL[NNsta][NNsta]; /*Left matrix quaternion*/

    double P[NEsta][NEsta];  /* prediction error covariance */
    double Q[NEsta][NEsta];  /* process noise covariance */
    double R[Mobs][Mobs];  /* measurement error covariance */

    double K[NEsta][Mobs];  /* Kalman gain; a.k.a. K */
    double Kt[Mobs][NEsta];  /* transpose Kalman gain; a.k.a. K */

    double Fx[NNsta][NNsta];  /* Jacobian of process model */
    double Fdx[NEsta][NEsta];  /* Jacobian of process model */
    double H[Mobs][NEsta];  /* Jacobian of measurement model */

    double Ht[NEsta][Mobs]; /* transpose of measurement Jacobian */
    double Fdxt[NEsta][NEsta]; /* transpose of process Jacobian */
    double Pp[NEsta][NEsta]; /* P, post-prediction, pre-update */
    
    double G[NEsta][NEsta];  

    double fx[NNsta];   /* output of user defined f() state-transition function */
    double hx[Mobs];   /* output of user defined h() measurement function */

    /* temporary storage */
    double tmp0[NEsta][NEsta];
    double tmp1[NEsta][Mobs];
    double tmp2[Mobs][NEsta];
    double tmp3[Mobs][Mobs];
    double tmp4[Mobs][Mobs];
    double tmp5[Mobs];
    double tmp6[NNsta]; 
    double tmp7[NNsta];

} ekf_t;        
