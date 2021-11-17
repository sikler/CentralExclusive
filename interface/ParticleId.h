#ifndef _ParticleId_h_
#define _ParticleId_h_

#include <vector>

#include "../interface/Parameters.h"
#include "../interface/Histo.h"

struct CmTrack;
struct Vector3; 

class MostProbable;
class Random;

class ParticleId
{
 public:
  ParticleId(int flag);
  ~ParticleId();

  //
  void processEnergyLoss(const std::vector<CmTrack> & tracks,
                         int evCategory, int type);

  // return particle type (unknown, pion, kaon, proton) for data
  int  identifyPair(const std::vector<CmTrack> & tra);

  // return type // // isIdentified for simulation
  int identifyPair(const std::vector<Vector3> & ps, int type);

  void prepareDemo(int type);

  double getEfficiency(double eta3, double p3t,
                       double eta4, double p4t, int type);

 private:
  double sampleRelSigma(double eta, double pt, const Histo & his);

  void propeller(int i);

  MostProbable * mostProbable;
  Random       * theRandom; 

  Histo elossAl,elossNo,elossPi,elossKa,elossPr, elossSb;
  Histo relSigma[nParts];
};

#endif
