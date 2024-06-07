/**************************************************************************
* Authors: 
*  Matt Benson and Daniel Malasky
* Summary:
*  Physics for the Artillery Prototype
**************************************************************************/

#define _USE_MATH_DEFINES
#include <math.h>   // for M_PI which is 3.14159

//#include "position.h"
#include <utility>
#include <vector>
#include <cmath>    // for sin, cos
#include <iostream>
using namespace std;

#define WEIGHT    46.7   // weight of projectile
#define DIAMETER  154.89 // diameter of projectile
#define VELOCITYM 827.0  // muzzle velocity
#define TIME      0.01   // time interval
#define GRAVITY   -9.8   // gravity
#define ACC       43.867 // acceleration


/****************************************************
* Table values
****************************************************/
// Altitude values
std::vector<int> altitudes =
{
 0, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000,
 10000, 15000, 20000, 25000, 30000, 40000, 50000, 60000, 70000, 80000
};

// Speed of Sound Values
std::vector<int> speedOfSound =
{
 340, 336, 332, 328, 324, 320, 316, 312, 308, 303,
 299, 295, 295, 295, 305, 324, 337, 319, 289, 269
};

// Density Values
std::vector<double> densities =
{
 1.2250000, 1.1120000, 1.0070000, 0.9093000, 0.8194000, 0.7364000,
 0.6601000, 0.5900000, 0.5258000, 0.4671000, 0.4135000, 0.1948000,
 0.0889100, 0.0400800, 0.0184100, 0.0039960, 0.0010270, 0.0003097,
 0.0000828, 0.0000185
};

// Mach values
std::vector<double> machNumbers =
{
 0.300, 0.500, 0.700, 0.890, 0.920, 0.960, 0.980, 1.000, 1.020, 1.060,
 1.240, 1.530, 1.990, 2.870, 2.890, 5.000
};

// Drag Coefficients
std::vector<double> dragCoefficients =
{
 0.1629, 0.1659, 0.2031, 0.2597, 0.3010, 0.3287, 0.4002, 0.4258, 0.4335, 0.4483,
 0.4064, 0.3663, 0.2897, 0.2297, 0.2306, 0.2656
};

/****************************************************
* Linear Interpolation
* Find a point between two points 
* d0,r0 = coordinates of one point
* d1,r1 = coordinates of second point
* d,r   = coordinates of a point in the middle
****************************************************/
double linearInterpolation(double d0, double r0, double d1, double r1, /*double d,*/ double r)
{
    
   return d0 + (((r - r0) * (d1 - d0)) / (r1 - r0));
}

/****************************************************
* ConvertToRadians
* Take an angle and convert it to radians
****************************************************/
double convertToRadians(double degrees)
{
   return degrees / 360 * (2 * M_PI);
}

/***************************************************
* COMPUTE DISTANCE
* Apply inertia to compute a new position using the distance equation.
* The equation is:
*     s = s + v t + 1/2 a t^2
* INPUT
*     s : original position, in meters
*     v : velocity, in meters/second
*     a : acceleration, in meters/second^2
*     t : time, in seconds
* OUTPUT
*     s : new position, in meters
**************************************************/
double computeDistance(double initialS, double v, double a, double t)
{
   double s;
   s = initialS + (v * t) + (0.5 * a * t * t);
   return s;
}

/****************************************************
* ComputeDX
****************************************************/
double computeDX(double angle, double initialV)
{
   return sin(angle) * initialV;
}

/****************************************************
* ComputeDY
****************************************************/
double computeDY(double angle, double initialV)
{
   return cos(angle) * initialV;
}


/****************************************************
* ComputeDDX
****************************************************/
double computeDDX(double angle)
{
   return (-sin(angle) * ACC);
}

/****************************************************
* ComputeDDY
****************************************************/
double computeDDY(double angle)
{
   return (GRAVITY - cos(angle) * ACC);
}


int main()
{
   pair<double, double> position(0, 0);
   double degrees = 75;
   double radians = convertToRadians(75);
   double radius = DIAMETER / 2;
   double surfaceA = M_PI * (radius * radius);

   //cout << "What is the angle of the howitzer where 0 is up?";
   //cin >> angle;

   cout << (linearInterpolation(0.2897, 1.900, 0.2297, 2.870, 2.4323));

   double distance = 0;
   double altitude = 0;
   double dx, dy;

   for (int i = 0; i < 20; i++)
   {
      double angle = convertToRadians(75);
      
      dx = computeDX(angle, VELOCITYM);
      dy = computeDY(angle, VELOCITYM);

      double ddx = computeDDX(angle);
      double ddy = computeDDY(angle);

      distance = computeDistance(0.0, dx, ddx, 20.0);
      altitude = computeDistance(0.0, dy, ddy, 20.0);
   }

   cout << distance << endl;
   cout << altitude << endl;


}


// Ask about inertia

// distance formula is inertia

// gravity - use linear interpolation and a return a value from a table
// drag - same but new table. drag force equatoin


// include position.h, velocity.h, and angle.h
