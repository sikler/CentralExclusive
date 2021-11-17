#ifndef _Physics_h_
#define _Physics_h_

//
#include <vector>

#include "../interface/Parameters.h"

//
struct Vector;
struct Vector3;
struct Vector4;
struct Event;

class RpVeto;
class TrkReco;

class Distributions;
class Resolutions;

typedef std::pair<double,Vector3> Boost; // FIXME

//
class Physics
{
 public:
  Physics(int nCmTra);
  ~Physics();

  void process(const Event & event);

 private:
  void correctMomenta(const Event & event,
                      Vector4 & pa, Vector4 & pb,
                      Vector4 & p1, Vector4 & p2, std::vector<Vector4> & hs);

  void getAngles(const Vector3 & ay, const Vector3 & az,
                 const Vector3 & p, std::vector<double> & a);

  void calculateAngles(const Vector4 & pa, const Vector4 & pb,
                             Vector4 q1, Vector4 q2,
                             Vector4 p3, Vector4 p4,
                             std::vector<double> & a);

  double getPhi_tilde(Vector4 p1, Vector4 p2, Vector4 q1,
                      const Boost & b_to_cm);


  void processMomenta(const Event & event,
                      const Vector4 & pa, const Vector4 & pb,
                      const Vector4 & p1, const Vector4 & p2,
                      const std::vector<Vector4> & hs);

  RpVeto  * theRpVeto;
  TrkReco * theTrkReco;

  Distributions * theDistributions;
  Resolutions   * theResolutions;

  double effIntLumi[nTopos];
};

#endif
