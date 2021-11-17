#ifndef _Distributions_h_
#define _Distributions_h_

#include <string>
#include <vector>

#include "../interface/Parameters.h"
#include "../interface/Histo.h"

#define nRegs 4 // FIXME

class Vector4;

class Distributions
{
 public:
  Distributions();
  ~Distributions();

  void collectForPhysics(int topo, int type,
       float dphi_tilde,
       const std::vector<float> & dphip1tp2t, float m,
       const Vector4 & q1, const Vector4 & q2,
       const Vector4 & p3, const Vector4 & p4, float weight);

 private:
  void average(   Histo his[], const std::vector<int> & ix, Histo  & ave);

  void symmetrize(Histo & his, const std::vector<int> & ix);

  // dphi - dphi_tilde
  Histo his_dphis[nParts];

  // regions
  double lowSign[nRegs],higSign[nRegs];
  int regType[nRegs];
  std::string regName[nRegs];

  // physics resonance
  Histo his_phi_topo[nRegs][nTopos];
  Histo his_phi[nRegs];

  Histo his_mph_topo[nParts][nTopos];
  Histo his_mph[nParts];

  Histo his_mas[nParts];
  Histo his_mhi[nParts];

  Histo his_sph_topo[nParts][nTopos];
  Histo his_sph[nParts];
  Histo his_sha[nParts]; // smaller of that and uhat

  // physics virtual
  Histo his_hat_topo[nParts][nTopos];
  Histo his_hat[nParts]; // that vs uhat
};

#endif
