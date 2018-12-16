/*
 * TinyEKF: Extended Kalman Filter for Arduino and TeensyBoard.
 *
 * Copyright (C) 2015 Simon D. Levy
 *
 * MIT License
 */

#include <stdio.h>
#include <stdlib.h>
#include "tiny_ekf_struct.h"

// Support both Arduino and command-line versions
#ifndef MAIN
extern "C" {
#endif
    void ekf_init(void *, int, int, int);
    int ekf_estimation(void *);
    int ekf_correction(void *, double *);
#ifndef MAIN
}
#endif

/**
 * A header-only class for the Extended Kalman Filter.  Your implementing class should #define the constant N and 
 * and then #include <TinyEKF.h>  You will also need to implement a model() method for your application.
 */
class TinyEKF {

    private:

        ekf_t ekf;

    protected:

        /**
          * The current state.
          */
        double * x;
        double t_lastCall;
        double dt;

        /**
         * Initializes a TinyEKF object.
         */
        TinyEKF() { 
            ekf_init(&this->ekf, NNsta, NEsta, Mobs);
            this->x = this->ekf.x; 
        }

        /**
         * Deallocates memory for a TinyEKF object.
         */
        ~TinyEKF() { }

        /**
         * Implement this function for your EKF model.
         * @param fx gets output of state-transition function <i>f(x<sub>0 .. n-1</sub>)</i>
         * @param F gets <i>n &times; n</i> Jacobian of <i>f(x)</i>
         * @param hx gets output of observation function <i>h(x<sub>0 .. n-1</sub>)</i>
         * @param H gets <i>m &times; n</i> Jacobian of <i>h(x)</i>
         */         
        virtual void model_estimation(double Fx[NNsta][NNsta], double Fdx[NEsta][NEsta], double meas[NEsta]) = 0;
        
        virtual void model_correction(double H[Mobs][NEsta], double x[NNsta], double h[Mobs], double qL[NNsta][NNsta]) = 0;
      
        virtual void setFx(double Fx[NNsta][NNsta], double meas[Mobs], double dt) = 0;
         
        virtual void setFdx(double Fdx[NEsta][NEsta], double meas[Mobs], double dt) = 0;
        
        virtual void setH(double H[NEsta][NEsta], double x[NNsta]) = 0;
        
        virtual void seth(double h[Mobs], double x[NNsta]) = 0;
        
        virtual void setqL(double qL[NNsta][NNsta], double x[NNsta]) = 0;
        /**
         * Sets the specified value of the prediction error covariance. <i>P<sub>i,j</sub> = value</i>
         * @param i row index
         * @param j column index
         * @param value value to set
         */
        void setP(int i, int j, double value) 
        { 
            this->ekf.P[i][j] = value; 
        }

        /**
         * Sets the specified value of the process noise covariance. <i>Q<sub>i,j</sub> = value</i>
         * @param i row index
         * @param j column index
         * @param value value to set
         */
        void setQ(int i, int j, double value) 
        { 
            this->ekf.Q[i][j] = value; 
        }

        /**
         * Sets the specified value of the observation noise covariance. <i>R<sub>i,j</sub> = value</i>
         * @param i row index
         * @param j column index
         * @param value value to set
         */
        void setR(int i, int j, double value) 
        { 
            this->ekf.R[i][j] = value; 
        }

    public:

        /**
         * Returns the state element at a given index.
         * @param i the index (at least 0 and less than <i>n</i>
         * @return state value at index
         */
        double getX(int i) 
        { 
            return this->ekf.x[i]; 
        }

        /**
         * Sets the state element at a given index.
         * @param i the index (at least 0 and less than <i>n</i>
         * @param value value to set
         */
        void setX(int i, double value) 
        { 
            this->ekf.x[i] = value; 
        }

        /**
          Performs one step of the prediction and update.
         * @param z observation vector, length <i>m</i>
         * @return true on success, false on failure caused by non-positive-definite matrix.
         */
        bool step(double * z) 
        { 
            double z_estim[3];
            double z_corre[3];
            
            z_estim[0] = z[0];
            z_estim[1] = z[1];
            z_estim[2] = z[2];
            
            z_corre[0] = z[3];
            z_corre[1] = z[4];
            z_corre[2] = z[5];
            
            this->model_estimation(this->ekf.Fx, this->ekf.Fdx, z_estim); 
            ekf_estimation(&this->ekf);

            // No sé segur si això funcionarà. Faltarà que el valor ekf.fx s'actualitzi des de "tiny_ekf.c"
            //this->model_correction(this->ekf.H, this->ekf.fx, this->ekf.hx, this->ekf.qL);
                      
            //return ekf_correction(&this->ekf, z_corre) ? false : true;
            return true;
        }
};
