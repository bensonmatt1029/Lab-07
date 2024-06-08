/**************************************************************************
* Authors: 
*  Matt Benson and Daniel Malasky
* Summary:
*  Physics for the Artillery Prototype
**************************************************************************/

#define _USE_MATH_DEFINES
#include <math.h>   // for M_PI which is 3.14159


#include "velocity.h"
#include "acceleration.h"
#include "position.h"
#include <utility>
#include <vector>
#include <cmath>    // for sin, cos
#include <iostream>
using namespace std;

#define WEIGHT    46.7   // weight of projectile
#define DIAMETER  154.89 // diameter of projectile
#define VELOCITYM 827.0  // muzzle velocity
#define TIME      0.01   // time interval
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
 0.1629, 0.1659, 0.2031, 0.2597, 0.3010, 0.3287, 0.4002, 0.4258, 0.4335,
 0.4483, 0.4064, 0.3663, 0.2897, 0.2297, 0.2306, 0.2656
};

// Gravity values
std::vector<double> gravities = 
{
 9.807, 9.804, 9.801, 9.797, 9.794, 9.791,
 9.788, 9.785, 9.782, 9.779, 9.776, 9.761, 9.745, 9.730
};


/****************************************************
* Linear Interpolation
* Find a point between two points 
* x0,y0 = coordinates of one point
* x1,y1 = coordinates of second point
* x,y   = coordinates of a point in the middle
****************************************************/
double linearInterpolation(double x0, double y0, double x1, double y1, double x) 
{
   return y0 + (y1 - y0) * ((x - x0) / (x1 - x0));
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

///******************************************************
//* setGravity
//* uses interpolation and values from the table
//******************************************************/
//int getIAltitudes(double altitude, int iAltitudes)
//{
//   
//   
//   //cout << "ALEVEL: " << aLevel << endl;
//   if (altitude > altitudes[iAltitudes+1] )
//   {
//      iAltitudes++;
//      
//   }
//   else if (altitude < altitudes[iAltitudes])
//   {
//      iAltitudes--;
//     
//   }
//   return iAltitudes;
//}

/******************************************************
* Get index of altitude in the altitude table
* Returns the index of the closest altitude to the given altitude
******************************************************/
int getAltitudeIndex(double altitude)
{
   int index = 0;
   while (index < altitudes.size() - 1 && altitude > altitudes[index + 1]) {
      index++;
   }
   return index;
}

/******************************************************
* Compute gravity at a given altitude
* Uses linear interpolation to find gravity at the given altitude
******************************************************/
double computeGravity(double altitude)
{
   int index = getAltitudeIndex(altitude);
   if (index == altitudes.size() - 1) {
      return gravities[index];
   }
   else {
      double gravity = linearInterpolation(altitudes[index], gravities[index], altitudes[index + 1], gravities[index + 1], altitude);
      return gravity;
   }
}

void setGravity(double altitude, int iAltitude)
{

}

int main()
{
   Acceleration a;

   double angle = convertToRadians(75);
   Velocity velocity(sin(angle) * VELOCITYM, cos(angle) * VELOCITYM);

   Position position;
   double degrees = 75;
   double radians = convertToRadians(75);
   double radius = DIAMETER / 2.0;
   double surfaceA = M_PI * (radius * radius);

   //cout << "What is the angle of the howitzer where 0 is up?";
   //cin >> angle;

   //cout << (linearInterpolation(0.2897, 1.900, 0.2297, 2.870, 2.4323));

   double distance = 0.0;
   double altitude = 0.0;
   double hangTime = 0.0;


   double gravity = -9.807;
   
   // leves on the grav table
   int iAltitudes = 0;
   int iGravities = 0;

   //for (int i = 0; i < 20; i++)
   while (altitude >= 0.0)
   {

      //iAltitudes = getIAltitudes(altitude, iAltitudes);
      //iGravities = iAltitudes;
      

      gravity = computeGravity(altitude);
      cout << "GRAV: " << gravity << endl;
      

      velocity.addDY(-gravity * TIME);

      distance += computeDistance(0.0, velocity.getDX(), 0.0, TIME);
      altitude += computeDistance(0.0, velocity.getDY(), -gravity, TIME);

      cout << "Distance: " << distance << " meters\n";
      cout << "Altitude: " << altitude << " meters\n";
      cout << "Velocity DY: " << velocity.getDY() << " m/s\n\n";
      cout << "Hang Time: " << hangTime << "s\n\n";

      hangTime += 0.01;
      
   }




}


// Ask about inertia

// distance formula is inertia

// gravity - use linear interpolation and a return a value from a table
// drag - same but new table. drag force equatoin


// include position.h, velocity.h, and angle.h
