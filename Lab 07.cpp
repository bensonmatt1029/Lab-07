/**************************************************************************
* Authors: 
*  Matt Benson and Daniel Malasky
* Summary:
*  Physics for the Artillery Prototype
**************************************************************************/

#define _USE_MATH_DEFINES
#include <math.h>             // for M_PI 
#include "velocity.h"         // for velocity
#include <vector>             // for vector
#include <cmath>              // for sin, cos
#include <iostream>           // for iostream
using namespace std;

#define WEIGHT    46.7        // weight of projectile
#define VELOCITYM 827.0       // muzzle velocity
#define TIME      0.01        // time interval
#define SURFACEAREA  0.018842 //in m^2 for drag

/****************************************************
* Table values
* Altitude, Speed of Sound, Air Density, Mach values,
* Drag Coefficients, Gravity Values
****************************************************/
// Altitude values
std::vector<double> altitudes =
{
 0, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000,
 10000, 15000, 20000, 25000, 30000, 40000, 50000, 60000, 70000, 80000
};

// Speed of Sound Values
std::vector<double> speedOfSound =
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

/*********************************************************************
* COMPUTE DRAG
* drag coefficient, airDensity, magnitude of Velocity, surface area
*********************************************************************/
double computeDrag(double c, double p, double v, double a)
{
   return (0.5 * c * p * (v * v) * a);
}

/****************************************************
* Linear Interpolation
* Find a point between two points 
* x0,y0 = coordinates of one point
* x1,y1 = coordinates of second point
* x,y   = coordinates of a point in the middle
****************************************************/
double linearInterpolation(double x0, double y0,
   double x1, double y1, double x) 
{
   return y0 + (y1 - y0) * ((x - x0) / (x1 - x0));
}

/****************************************************
* CONVERT TO RADAINS
* Take an angle and convert it to radians
****************************************************/
double convertToRadians(double degrees)
{
   return degrees / 360 * (2 * M_PI);
}

/********************************************************************
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
*******************************************************************/
double computeDistance(double initialS, double v, double a, double t)
{
   double s;
   s = initialS + (v * t) + (0.5 * a * t * t);
   return s;
}

/******************************************************
* GET TABLE INDEX
* Get index of a given value in the given table
* Returns the index of the closest values to the given value
******************************************************/
int getTableIndex(double value, vector<double> values)
{
   int index = 0;
   while (index < values.size() - 1 && value > values[index + 1]) 
   {
      index++;
   }
   return index;
}

/****************************************************************
* COMPUTE GRAVITY
* Compute gravity at a given altitude
* Uses linear interpolation to find gravity at the given altitude
*****************************************************************/
double computeGravity(double altitude)
{
   int index = getTableIndex(altitude, altitudes);
   if (index == altitudes.size() - 1) 
   {
      return gravities[index];
   }
   else 
   {
      double gravity = linearInterpolation(altitudes[index],
      gravities[index], altitudes[index + 1], gravities[index + 1], altitude);
      return gravity;
   }
}

/****************************************************************
* COMPUTE SPEED OF SOUND
* Uses linear interpolation to find speed of sound at the given altitude
*****************************************************************/
double computeSpeedOfSound(double altitude)
{
   int index = getTableIndex(altitude, altitudes);
   if (index == altitudes.size() - 1) 
   {
      return speedOfSound[index];
   }
   else 
   {
      double soundSpeed = linearInterpolation(altitudes[index], 
         speedOfSound[index], altitudes[index + 1], 
         speedOfSound[index + 1], altitude);
      return soundSpeed;
   }
}

/******************************************************
* COMPUTE DRAG COEFFICIENT
* Uses linear interpolation to find drag at the given mach
******************************************************/
double computeDragCo(double mach)
{
   int index = getTableIndex(mach, machNumbers);
   if (index == machNumbers.size() - 1) 
   {
      return machNumbers[index];
   }
   else 
   {
      double dragCo = linearInterpolation(machNumbers[index], 
         dragCoefficients[index], machNumbers[index + 1], 
         dragCoefficients[index + 1], mach);
      return dragCo;
   }
}

/***************************************************************
* COMPUTE AIR DENSITY
* Uses linear interpolation to find gravity at the given altitude
****************************************************************/
double computeAirDensity(double altitude)
{
   int index = getTableIndex(altitude, altitudes);
   if (index == densities.size() - 1) 
   {
      return densities[index];
   }
   else 
   {
      double airDensity = linearInterpolation(altitudes[index], 
      densities[index],altitudes[index + 1], densities[index + 1], altitude);
      return airDensity;
   }
}

/*************************************************************
* DISPLAY BULLET TRAVEL
* Displays the distance and hang time of a 
* M795 High Explosive 155mm artillery projectile
* from the M777 howitzer
******************************************************/
void displayBulletTravel(double angleDegrees) 
{
   double angle = convertToRadians(angleDegrees);
   Velocity velocity(sin(angle) * VELOCITYM, cos(angle) * VELOCITYM);

   double distance = 0.0;
   double altitude = 0.0;
   double hangTime = 0.0;
   double gravity;
   double airDensity;
   double dragCo;
   double soundSpeed;

   // second to last values
   double distance0 = 0;
   double altitude0 = 0;

   while (altitude >= 0.0) {
      // Compute values 
      gravity = computeGravity(altitude);
      airDensity = computeAirDensity(altitude);
      soundSpeed = computeSpeedOfSound(altitude);
      dragCo = computeDragCo(velocity.getSpeed() / soundSpeed);

      // Calculate drag force
      double totalVelocity = velocity.getSpeed();
      double dragForce = computeDrag(dragCo, airDensity,
         totalVelocity, SURFACEAREA);
      double dragAcceleration = dragForce / WEIGHT;

      // Calculate components of drag acceleration
      double dragAccX = -dragAcceleration * (velocity.getDX() / totalVelocity);
      double dragAccY = -dragAcceleration * (velocity.getDY() / totalVelocity);

      // Update velocities with drag and gravity
      velocity.addDX(dragAccX * TIME);
      velocity.addDY((-gravity + dragAccY) * TIME);

      // Update positions
      distance += computeDistance(0.0, velocity.getDX(), 0.0, TIME);
      altitude += computeDistance(0.0, velocity.getDY(), -gravity, TIME);

      // update hangTime
      hangTime += TIME;

      
      // second to last distance, altitude for Ground Impact
      if (altitude > 0)
      {
         distance0 = distance;
         altitude0 = altitude;
      }
   }
   
   // Ground Impact
   distance = linearInterpolation(altitude0, distance0, altitude, distance, 0);

   cout << "Distance: " << distance << " meters\t";
   cout << "Hang Time: " << hangTime << "s\n\n";
}

/*************************************************************
* MAIN
* There can be only one!
******************************************************/
int main()
{
   // Test Cases
   displayBulletTravel(0);
   displayBulletTravel(30);
   displayBulletTravel(60);
   displayBulletTravel(-45);
}