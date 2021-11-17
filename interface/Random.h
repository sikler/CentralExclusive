#ifndef _Random_h_
#define _Random_h_

class Random
{
 public:
  Random() {};
  ~Random() {};

  double getFlat();
  double getFlat(double a, double b);
  double getFlat(double a, double b, double r);

  double getGauss();
  double getLog();

 private:
};

#endif
