/*
 * TinyEKF: Extended Kalman Filter for embedded processors
 *
 * Copyright (C) 2015 Simon D. Levy
 *
 * MIT License
 */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

/* Cholesky-decomposition matrix-inversion code, adapated from
   http://jean-pierre.moreau.pagesperso-orange.fr/Cplus/choles_cpp.txt */

static int choldc1(double * a, double * p, int n) {
    int i,j,k;
    double sum;

    for (i = 0; i < n; i++) {
        for (j = i; j < n; j++) {
            sum = a[i*n+j];
            for (k = i - 1; k >= 0; k--) {
                sum -= a[i*n+k] * a[j*n+k];
            }
            if (i == j) {
                if (sum <= 0) {
                    return 1; /* error */
                }
                p[i] = sqrt(sum);
            }
            else {
                a[j*n+i] = sum / p[i];
            }
        }
    }

    return 0; /* success */
}

static int choldcsl(double * A, double * a, double * p, int n) 
{
    int i,j,k; double sum;
    for (i = 0; i < n; i++) 
        for (j = 0; j < n; j++) 
            a[i*n+j] = A[i*n+j];
    if (choldc1(a, p, n)) return 1;
    for (i = 0; i < n; i++) {
        a[i*n+i] = 1 / p[i];
        for (j = i + 1; j < n; j++) {
            sum = 0;
            for (k = i; k < j; k++) {
                sum -= a[j*n+k] * a[k*n+i];
            }
            a[j*n+i] = sum / p[j];
        }
    }

    return 0; /* success */
}


static int cholsl(double * A, double * a, double * p, int n) 
{
    int i,j,k;
    if (choldcsl(A,a,p,n)) return 1;
    for (i = 0; i < n; i++) {
        for (j = i + 1; j < n; j++) {
            a[i*n+j] = 0.0;
        }
    }
    for (i = 0; i < n; i++) {
        a[i*n+i] *= a[i*n+i];
        for (k = i + 1; k < n; k++) {
            a[i*n+i] += a[k*n+i] * a[k*n+i];
        }
        for (j = i + 1; j < n; j++) {
            for (k = j; k < n; k++) {
                a[i*n+j] += a[k*n+i] * a[k*n+j];
            }
        }
    }
    for (i = 0; i < n; i++) {
        for (j = 0; j < i; j++) {
            a[i*n+j] = a[j*n+i];
        }
    }

    return 0; /* success */
}

static void zeros(double * a, int m, int n)
{
    int j;
    for (j=0; j<m*n; ++j)
        a[j] = 0;
}

#ifdef DEBUG
static void dump(double * a, int m, int n, const char * fmt)
{
    int i,j;

    char f[100];
    sprintf(f, "%s ", fmt);
    for(i=0; i<m; ++i) {
        for(j=0; j<n; ++j)
            printf(f, a[i*n+j]);
        printf("\n");
    }
}
#endif

/* C <- A * B */
static void mulmat(double * a, double * b, double * c, int arows, int acols, int bcols)
{
    int i, j,l;

    for(i=0; i<arows; ++i)
        for(j=0; j<bcols; ++j) {
            c[i*bcols+j] = 0;
            for(l=0; l<acols; ++l)
                c[i*bcols+j] += a[i*acols+l] * b[l*bcols+j];
        }
}

static void mulvec(double * a, double * x, double * y, int m, int n)
{
    int i, j;

    for(i=0; i<m; ++i) {
        y[i] = 0;
        for(j=0; j<n; ++j)
            y[i] += x[j] * a[i*n+j];
    }
}

static void transpose(double * a, double * at, int m, int n)
{
    int i,j;

    for(i=0; i<m; ++i)
        for(j=0; j<n; ++j) {
            at[j*m+i] = a[i*n+j];
        }
}

/* A <- A + B */
static void accum(double * a, double * b, int m, int n)
{        
    int i,j;

    for(i=0; i<m; ++i)
        for(j=0; j<n; ++j)
            a[i*n+j] += b[i*n+j];
}

/* C <- A + B */
static void add(double * a, double * b, double * c, int n)
{
    int j;

    for(j=0; j<n; ++j)
        c[j] = a[j] + b[j];
}


/* C <- A - B */
static void sub(double * a, double * b, double * c, int n)
{
    int j;

    for(j=0; j<n; ++j)
        c[j] = a[j] - b[j];
}

static void negate(double * a, int m, int n)
{        
    int i, j;

    for(i=0; i<m; ++i)
        for(j=0; j<n; ++j)
            a[i*n+j] = -a[i*n+j];
}

static void mat_addeye(double * a, int n)
{
    int i;
    for (i=0; i<n; ++i)
        a[i*n+i] += 1;
}

/* TinyEKF code ------------------------------------------------------------------- */

#include "tiny_ekf.h"

typedef struct {

    double * x;    /* state vector */
    double * dx;   /* error-state vector*/

    double * P;  /* prediction error covariance */
    double * Q;  /* process noise covariance */
    double * R;  /* measurement error covariance */

    double * G;  /* Kalman gain; a.k.a. K */

    double * Fx;  /* Jacobian of process model */
    double * Fdx;  /* Jacobian of process model */
    double * H;  /* Jacobian of measurement model */

    double * Ht; /* transpose of measurement Jacobian */
    double * Fdxt; /* transpose of process Jacobian */
    double * Pp; /* P, post-prediction, pre-update */

    double * fx;  /* output of user defined f() state-transition function */
    double * hx;  /* output of user defined h() measurement function */

    /* temporary storage */
    double * tmp0;
    double * tmp1;
    double * tmp2;
    double * tmp3;
    double * tmp4;
    double * tmp5; 

} ekf_t;

static void unpack(void * v, ekf_t * ekf, int nn, int ne, int m)
{
    /* skip over nn, ne, m in data structure */
    char * cptr = (char *)v;
    cptr += 3*sizeof(int);

    double * dptr = (double *)cptr;
    ekf->x = dptr;
    dptr += nn;
    ekf->dx = dptr;
    dptr += ne;
    ekf->P = dptr;
    dptr += ne*ne;
    ekf->Q = dptr;
    dptr += ne*ne;
    ekf->R = dptr;
    dptr += m*m;
    ekf->G = dptr;
    dptr += ne*m;
    ekf->Fx = dptr;
    dptr += nn*nn;
    ekf->Fdx = dptr;
    dptr += ne*ne;
    ekf->H = dptr;
    dptr += m*ne;
    ekf->Ht = dptr;
    dptr += ne*m;
    ekf->Fdxt = dptr;
    dptr += ne*ne;
    ekf->Pp = dptr;
    dptr += ne*ne;
    ekf->fx = dptr;
    dptr += nn;
    ekf->hx = dptr;
    dptr += m;
    ekf->tmp0 = dptr;
    dptr += ne*ne;
    ekf->tmp1 = dptr;
    dptr += ne*m;
    ekf->tmp2 = dptr;
    dptr += m*ne;
    ekf->tmp3 = dptr;
    dptr += m*m;
    ekf->tmp4 = dptr;
    dptr += m*m;
    ekf->tmp5 = dptr;
  }

void ekf_init(void * v, int nn, int ne, int m)
{
    /* retrieve nn, ne, m and set them in incoming data structure */
    /* Set value to the ekf object passed as v*/
    int * ptr = (int *)v;
    *ptr = nn;
    ptr++;
    *ptr = ne;
    ptr++;
    *ptr = m;

    /* unpack rest of incoming structure for initlization */
    ekf_t ekf;
    unpack(v, &ekf, nn, ne, m);

    /* zero-out matrices */
    zeros(ekf.P, ne, ne);
    zeros(ekf.Q, ne, ne);
    zeros(ekf.R, m, m);
    zeros(ekf.G, ne, m);
    zeros(ekf.Fx, nn, nn);
    zeros(ekf.Fdx, ne, ne);
    zeros(ekf.H, m, ne);
}

int ekf_step(void * v, double * z)
{        
    /* unpack incoming structure */

    int * ptr = (int *)v;
    int nn = *ptr;
    ptr++;
    int ne = *ptr;
    ptr++;
    int m = *ptr;

    ekf_t ekf;
    unpack(v, &ekf, nn, ne, m); 
    
    /* Remember to predict here the new state and store it in f(x), as it was 
    done before in the model method /*
    /* f(x) = F*ekf.x; */
    mulvec(ekf.Fx, ekf.x, ekf.fx, nn, nn);
    
    /* P_k = Fdx_{k-1} P_{k-1} Fdx^T_{k-1} + Q_{k-1} */
    mulmat(ekf.Fdx, ekf.P, ekf.tmp0, ne, ne, ne);
    transpose(ekf.Fdx, ekf.Fdxt, ne, ne);
    mulmat(ekf.tmp0, ekf.Fdxt, ekf.Pp, ne, ne, ne);
    accum(ekf.Pp, ekf.Q, ne, ne);
    
    /* Compute the residual */
    
    /* G_k = P_k H^T_k (H_k P_k H^T_k + R)^{-1} */
    transpose(ekf.H, ekf.Ht, m, ne);
    mulmat(ekf.Pp, ekf.Ht, ekf.tmp1, ne, ne, m);
    mulmat(ekf.H, ekf.Pp, ekf.tmp2, m, ne, ne);
    mulmat(ekf.tmp2, ekf.Ht, ekf.tmp3, m, ne, m);
    accum(ekf.tmp3, ekf.R, m, m);
    if (cholsl(ekf.tmp3, ekf.tmp4, ekf.tmp5, m)) return 1;
    mulmat(ekf.tmp1, ekf.tmp4, ekf.G, ne, m, m);

    /* \hat{x}_k = \hat{x_k} + G_k(z_k - h(\hat{x}_k)) */
    sub(z, ekf.hx, ekf.tmp5, m);
    mulvec(ekf.G, ekf.tmp5, ekf.tmp2, ne, m);
    add(ekf.fx, ekf.tmp2, ekf.x, ne);

    /* P_k = (I - G_k H_k) P_k */
    mulmat(ekf.G, ekf.H, ekf.tmp0, ne, m, ne);
    negate(ekf.tmp0, ne, ne);
    mat_addeye(ekf.tmp0, ne);
    mulmat(ekf.tmp0, ekf.Pp, ekf.P, ne, ne, ne);

    /* success */
    return 0;
}
