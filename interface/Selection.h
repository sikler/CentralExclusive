#ifndef _Selection_h_
#define _Selection_h_

//
#include "../interface/Parameters.h"

#include "../interface/Histo.h"

struct Event;
struct Vector2;

//
class Selection
{
 public:
  Selection(int nCmTra);
  ~Selection();

  bool     areCmTracksLooper(const Event & event, int type, double w);
  bool  areCmTracksPrimaries(const Event & event, double w);
  bool areRpProtonsPrimaries(const Event & event, double w);

  int classifyEvent(const Event & event, int topo, double w);

 private:
  double getChi(const Vector2 & v);

  void sumMomenta(const Event & event, Vector2 & sum2pt, Vector2 & sum4pt);

  Histo his_cm_xy, his_cm_zr;
  Histo his_cm_z, his_cm_dz, his_cm_dzrel, his_cm_dz_vs_sz;

  Histo his_cm_looper[nParts];

  Histo his_rp_x12_elas, his_rp_x12_sign, his_rp_x12_side,
        his_rp_dxrel, his_rp_dx_vs_sx;

  Histo his_chi_2, his_chi_4, his_chi_24;

  Histo his_chi_4_topo[nTopos];

  Histo his_oth_2, his_oth_4;

  Histo his_chi[nDphi];
  Histo his_oth[nDphi];
  Histo his_chi_24_dphi[nDphi];

  //
  Histo his_sum2p[nTopos], his_sum4p[nTopos],
        his_sumpx[nTopos], his_sumpy[nTopos];

  Histo his_sum2px_vs_pred[nTopos],
        his_sum2px_pred[nTopos],
        his_sum2py[nTopos];

  const int nCmTra;
};

#endif
