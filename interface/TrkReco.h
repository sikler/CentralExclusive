#ifndef _TrkReco_h_
#define _TrkReco_h_

#include <vector>
#include <map>

#include "../interface/Histo.h"
#include "../interface/HistoVal.h"
#include "../interface/Parameters.h"

//
const int minLayers   = 3;
const int minClusters = 5;

const int minBpixTracks  = 1;
const int minPixelTracks = 1;

// FIXME VERY
const int mPt = (5 * nPt)/2;

struct Part;
struct SimTrackInfo;

struct Vector3;
struct Vector4;
struct CmTrack;

class Resolutions;

//
class TrkReco
{
 public:
  TrkReco(int flag);
  ~TrkReco();

  void collectPars(const std::vector<CmTrack> & tracks);
  void collectChi2(const std::vector<CmTrack> & tracks, int type);

  void prepareCmTables(int type);

  bool hasGoodEffic(const std::vector<CmTrack> & tracks);

  void getMomenta(const Vector4 & pa, const Vector4 & pb,
                        Vector4   p1,       Vector4   p2,
                  int type, const std::vector<double> & ptmcosthetaphi,
                  double y,
                        Vector4 & p3,       Vector4 & p4);

  bool notFake(const Vector3 & p3,
               const Vector3 & p4, int type);

  double getCmEff(const std::vector<float> & dphip1tp2t,
                  const std::vector<float> & ptmcosthetaphi,
                  int topo, int type);

 private:
  bool hasGoodEffic(const Vector3 & p);

  void readCmEff(const int type);
  void readCmEffTables();

  HistoVal his_cm_eff[nTopos][nParts][nKt][nKt];
//  Histo his_cm_eff[nTopos][nParts][nKt][nKt];

  Histo his_trackZ, his_trackT, his_trackDz;
  Histo /*his_trackSp[nParts],*/ his_trackChi2[nParts], his_trackY[nParts];
  Histo his_resoY[nParts];

  Histo his_trackEtaPtPos,  his_trackEtaPtNeg;
  Histo his_trackEtaPhiPos, his_trackEtaPhiNeg;

  //
  Histo his_oneTrack_pos[nParts], his_oneTrack_neg[nParts];

  std::mutex mtx_cerr;
};

#endif
