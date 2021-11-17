#include "../interface/Random.h"

#include <cmath>

using namespace std;

/*****************************************************************************/
double Random::getFlat()
{
  return drand48();
}

/*****************************************************************************/
double Random::getFlat(double a, double b)
{
  return a + (b-a)*drand48();
}

/*****************************************************************************/
double Random::getFlat(double a, double b, double r)
{
  return a + (b-a)*r;
}

/*****************************************************************************/
double Random::getGauss()
{
  return sqrt(-2*log(getFlat())) * cos(M_PI*getFlat());
}

/*****************************************************************************/
double Random::getLog()
{
  return -log(getFlat());
}

