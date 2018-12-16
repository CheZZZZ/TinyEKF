/* SensorFusion: Sensor fusion on Arduino using TinyEKF.  
 *
 * Copyright (C) 2015 Simon D. Levy
 *
 * MIT License
 */


// These must be defined before including TinyEKF.h
#define Nsta 2     // Two state values: pressure, temperature
#define Mobs 3     // Three measurements: baro pressure, baro temperature, LM35 temperature

#define LM35_PIN 0

#include <TinyEKF.h>
#include <SFE_BMP180.h>
#include <Wire.h>

class Fuser : public TinyEKF {

    public:

        Fuser()/* SensorFusion: Sensor fusion on Arduino using TinyEKF.  
 *
 * Copyright (C) 2015 Simon D. Levy
 *
 * MIT License
 */


// These must be defined before including TinyEKF.h
#define NNsta 4     // Four nominal state values: quaternion
#define NEsta 3     // Three error state values: angle error
#define Mobs 3      // Six measurements: gyroscope, accelerometer. Only acc. couting as the correction

#include <Wire.h>
#include <math.h>
#include "LSM6DSM.h"
#include "TinyEKF.h"

//LSM6DSM definitions
#define INTERRUPT_PIN 2  // interrupt1 pin definitions, data ready
static LSM6DSM::Ascale_t Ascale = LSM6DSM::AFS_4G;
static LSM6DSM::Gscale_t Gscale = LSM6DSM::GFS_2000DPS;
static LSM6DSM::Rate_t AODR = LSM6DSM::ODR_1660Hz;
static LSM6DSM::Rate_t GODR = LSM6DSM::ODR_1660Hz;

// Gyro and accel parameters
static float GYRO_NOISE_DEN  = 3.8f;   // [mdps/sqrt(Hz)]
static float ACCEL_NOISE_DEN = 80.0f;  // [ug/sqrt(Hz)] (Different value depending on the FS)

// Random Bias to initiate the object
float ACCEL_BIAS[3] = {0.0f, 0.0f, 0.0f};
float GYRO_BIAS[3]  = {0.0f, 0.0f, 0.0f};

class Fuser : public TinyEKF {

    public:

        Fuser()
        {            
            // We approximate the process noise using a small constant
            this->setQ(0, 0, .0001);
            this->setQ(1, 1, .0001);

            // Same for measurement noise
            this->setR(0, 0, .0001);
            this->setR(1, 1, .0001);
            this->setR(2, 2, .0001);
            
        }

    protected:
              
        void model_estimation(double Fx[NNsta][NNsta], double Fdx[NEsta][NEsta], double meas[NEsta])
        {
            // Compute deltat
            double t_now = micros();
            this->dt = (t_now - this->t_lastCall)/1000000.0;
            this->t_lastCall = t_now;

            // Set nominal state Jacobian
            setFx(Fx, meas, this->dt);

            // Set error-sate Jacobian
            this->setFdx(Fdx, meas, this->dt);

        }

        void model_correction(double H[Mobs][NEsta], double x[NNsta], double h[Mobs], double qL[NNsta][NNsta])
        {
            // Set the jacobian for the correction measurement
            setH(H, x);

            // Set the measurement prediction for the estimated state
            seth(h,x);

            // Set the left quaternion matrix
            setqL(qL, x);
        }


        void setFx(double Fx[NNsta][NNsta], double meas[NEsta],double dt)
        {
            // First Column
            Fx[0][0]  = 1.0;
            Fx[1][0]  = meas[0]*dt/2.0;
            Fx[2][0]  = meas[1]*dt/2.0;
            Fx[3][0]  = meas[2]*dt/2.0;
            // Second Column
            Fx[0][1]  = -meas[0]*dt/2.0;
            Fx[1][1]  =  1.0;
            Fx[2][1]  = -meas[2]*dt/2.0;
            Fx[3][1]  =  meas[1]*dt/2.0;
            // Third Column
            Fx[0][2]  = -meas[1]*dt/2.0;
            Fx[1][2]  =  meas[2]*dt/2.0;
            Fx[2][2] =  1.0;
            Fx[3][2] = -meas[0]*dt/2.0;
            // Fourth Column
            Fx[0][3] = -meas[2]*dt/2.0;
            Fx[1][3] = -meas[1]*dt/2.0;
            Fx[2][3] =  meas[0]*dt/2.0;
            Fx[3][3] =  1.0;
        }

        void setFdx(double Fdx[NEsta][NEsta], double meas[NEsta],double dt)
        {
            // First Column
            Fdx[0][0] =  1.0;
            Fdx[1][0] =  meas[2]*dt;
            Fdx[2][0] = -meas[1]*dt;
            // Second Column
            Fdx[0][1] = -meas[2]*dt;
            Fdx[1][1] =  1.0;
            Fdx[2][1] =  meas[0]*dt;
            // Third Column
            Fdx[0][2] =  meas[1]*dt;
            Fdx[1][2] = -meas[0]*dt;
            Fdx[2][2] = 1.0;
        }

        void setH(double H[Mobs][NEsta], double x[NNsta])
        {
            // First Column
            H[0][0]  =  0.0;
            H[1][0]  =  x[0]*x[0] - x[1]*x[1] - x[2]*x[2] + x[3]*x[3];
            H[2][0]  = - 2*x[0]*x[1] - 2*x[2]*x[3];
            // Second Column
            H[0][1]  = -x[0]*x[0] + x[1]*x[1] + x[2]*x[2] - x[3]*x[3];
            H[1][1]  =  0.0;
            H[2][1]  = 2*x[1]*x[3] - 2*x[0]*x[2];
            // Third Column
            H[0][2]  = 2*x[0]*x[1] + 2*x[2]*x[3];
            H[1][2]  = 2*x[0]*x[2] - 2*x[1]*x[3];
            H[2][2]  = 0.0;
        }

        void seth(double h[Mobs], double x[NNsta])
        {
            // First Column
            h[0] = ( 2*x[0]*x[2] - 2*x[1]*x[3])*-9.80665;
            h[1] = (-2*x[0]*x[1] - 2*x[2]*x[3])*-9.80665;
            h[2] = (-x[0]*x[0] + x[1]*x[1] + x[2]*x[2] - x[3]*x[3])*-9.80665;
        }

        void setqL(double qL[NNsta][NNsta], double x[NNsta])
        {
            qL[0][0] =  x[0];
            qL[1][0] =  x[1];
            qL[2][0] =  x[2];
            qL[3][0] =  x[3];
            
            qL[0][1] = -x[1];
            qL[1][1] =  x[0];
            qL[2][1] =  x[3];
            qL[3][1] = -x[2];
            
            qL[0][2]  = -x[2];
            qL[1][2]  = -x[3];
            qL[2][2] =  x[0];
            qL[3][2] =  x[1];
            
            qL[0][3] = -x[3];
            qL[1][3] =  x[2];
            qL[2][3] = -x[1];
            qL[3][3] =  x[0];
        }

};

Fuser ekf;
static LSM6DSM lsm6dsm(Ascale, Gscale, AODR, GODR, ACCEL_BIAS, GYRO_BIAS);

void setup() {

    Serial.begin(115200);
   
    // Configure interrupt
    pinMode(INTERRUPT_PIN, INPUT);

    // Start I^2C 
    Wire.begin(TWI_PINS_20_21);           // set master mode 
    Wire.setClock(400000);  // I2C frequency at 400 kHz  
    delay(100);

    lsm6dsm.begin();
    
    IMUCalibrate();

    delay(1000);
}


void loop() {

    // Read IMU measurements
    if (lsm6dsm.checkNewData()) 
    {
      float ax=0.0f, ay=0.0f, az=0.0f, gx=0.0f, gy=0.0f, gz=0.0f;
      lsm6dsm.readData(ax, ay, az, gx, gy, gz);

      double z[6]  = {gx*M_PI/180.0, gy*M_PI/180.0, gz*M_PI/180.0, ax*9.80665, ay*9.80665, az*9.80665};           
      ekf.step(z);
    }
    
    // Report measured and predicte/fused values
    /*Serial.print(z[0]);
    Serial.print(" ");
    Serial.print(z[1]);
    Serial.print(" ");
    Serial.print(z[2]);
    Serial.print(" ");
    Serial.print(ekf.getX(0));
    Serial.print(" ");
    Serial.println(ekf.getX(1));*/
}

void IMUCalibrate()
{
  Serial.println("DO NOT move the IMU. Calibration starting in:");
  Serial.println("3");
  delay(500);
  Serial.println("2");
  delay(500);
  Serial.println("1");
  delay(500);
  
  Serial.println("Calibrating...");
  lsm6dsm.calibrate(GYRO_BIAS, ACCEL_BIAS);
  Serial.println("Calibrated");
}

        {            
            // We approximate the process noise using a small constant
            this->setQ(0, 0, .0001);
            this->setQ(1, 1, .0001);

            // Same for measurement noise
            this->setR(0, 0, .0001);
            this->setR(1, 1, .0001);
            this->setR(2, 2, .0001);
        }

    protected:

        void model(double fx[Nsta], double F[Nsta][Nsta], double hx[Mobs], double H[Mobs][Nsta])
        {
            // Process model is f(x) = x
            fx[0] = this->x[0];
            fx[1] = this->x[1];

            // So process model Jacobian is identity matrix
            F[0][0] = 1;
            F[1][1] = 1;

            // Measurement function simplifies the relationship between state and sensor readings for convenience.
            // A more realistic measurement function would distinguish between state value and measured value; e.g.:
            //   hx[0] = pow(this->x[0], 1.03);
            //   hx[1] = 1.005 * this->x[1];
            //   hx[2] = .9987 * this->x[1] + .001;
            hx[0] = this->x[0]; // Barometric pressure from previous state
            hx[1] = this->x[1]; // Baro temperature from previous state
            hx[2] = this->x[1]; // LM35 temperature from previous state

            // Jacobian of measurement function
            H[0][0] = 1;        // Barometric pressure from previous state
            H[1][1] = 1 ;       // Baro temperature from previous state
            H[2][1] = 1 ;       // LM35 temperature from previous state
        }
};

Fuser ekf;
SFE_BMP180 baro;

void setup() {

    Serial.begin(9600);

    // Start reading from baro
    baro.begin();

    // Set up to read from LM35
    analogReference(INTERNAL);
}


void loop() {

    // Read pressure, temperature from BMP180
    double baroTemperature, baroPressure;
    getBaroReadings(baroTemperature, baroPressure);

    // Read temperature from LM35
    float lm35Temperature = analogRead(LM35_PIN) / 9.31;

    // Send these measurements to the EKF
    double z[3] = {baroPressure, baroTemperature, lm35Temperature};
    ekf.step(z);

    // Report measured and predicte/fused values
    Serial.print(z[0]);
    Serial.print(" ");
    Serial.print(z[1]);
    Serial.print(" ");
    Serial.print(z[2]);
    Serial.print(" ");
    Serial.print(ekf.getX(0));
    Serial.print(" ");
    Serial.println(ekf.getX(1));
}


// Adapted from https://github.com/sparkfun/BMP180_Breakout
void getBaroReadings(double & T, double & P)
{
    char status = baro.startTemperature();
   
    if (status != 0) {
        delay(status);
        status = baro.getTemperature(T);
        if (status != 0) {
            status = baro.startPressure(3);
            if (status != 0) {
                delay(status);
                status = baro.getPressure(P,T);
                if (status == 0)
                    Serial.println("error retrieving pressure measurement");
            }
            else Serial.println("error starting pressure measurement");
        }
        else Serial.println("error retrieving temperature measurement");
    }
    else Serial.println("error starting temperature measurement");
}
