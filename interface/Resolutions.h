#ifndef _Resolutions_h_
#define _Resolutions_h_

#include "../interface/Parameters.h"

#include "../interface/Histo.h"

class Random;
struct Vector4;

class Resolutions
{
 public:
  Resolutions(int flag, const std::string & partName); // for simulation
  Resolutions(int flag);                               // for data

  ~Resolutions();

  void collectForPt(const std::vector<float> & etaptdpt);
  void collectForMass(int type, const std::vector<float> & dphip1tp2t,
                      double m, double weight,
                      const Vector4 & p3, const Vector4 & p4);

 private:
  Histo ptResol; // collect from simu

  Histo ptResol_mean[nParts][2], // for data correction, [charge]
        ptResol_sigm[nParts][2];

  Histo his_massResp[nKt][nKt][nParts];

  Random * theRandom;

  const int flag;
};

#endif
