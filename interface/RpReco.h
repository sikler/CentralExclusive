#ifndef _RpReco_h_
#define _RpReco_h_

#include <vector>
#include <map>
#include <mutex>

#include "../interface/Parameters.h"
#include "../interface/Histo.h"

//
class RpFit;
class RpEffic;

struct Vector2;
struct RpRes;
struct RpDet;

struct Event;
struct RpTrack;
struct PrTrack;

//
class RpReco
{
 public:
  RpReco(int nCmTra);
  ~RpReco();

  bool getTopo(const std::vector<RpTrack> & rpTracks, int & topo);

  bool reconstruct(int topo, Event & event); // isWithin

 private:
  void covariance(Histo & his,
                  double & sig_xstar, double & sig_tstar,
                  double & Ln, double & Lf);

  void localReco(int run, int topo, RpTrack & track,
                 std::map<RpDet,double> & locCy,
                 std::map<RpDet,double> & locCx);

  void correctPos(RpTrack & track);

  void getStar(int arm, const Vector2 pos[], RpRes & res);

  void globalReco(int topo,
                  const std::vector<RpTrack> & rpTracks,
                        std::vector<PrTrack> & prTracks,
                  double & rpWeight);

  RpFit   * theRpFit;
  RpEffic * theRpEffic;

  // optics
  double v_x[2],v_y[2], L_x[2][2],L_y[2];

  // shifts
  double mx[2][2][2], my[2][2][2];

  // histograms
  Histo his_weight[nTopos];

  Histo his_prot_loc[2][nTopos],
        his_prot_mom[2][nTopos],
        his_prot_mom_x[nTopos],
        his_prot_mom_y[nTopos];

  Histo his_vtx_x[nTopos], his_vtx_y[nTopos];

  //
  Histo his_hits_loc[2][2][2][nTopos],
        his_hits_loc_x[2][2][2][nTopos],
        his_hits_loc_y[2][2][2][nTopos];

  Histo his_hits_nf_x[2][nTopos],
        his_hits_nf_y[2][nTopos];

  //
  std::map<RpDet,Histo> his_occup;

  //
  int nCmTra;

  std::mutex mtx_rpfit, mtx_his, mtx_cout;
};

#endif
