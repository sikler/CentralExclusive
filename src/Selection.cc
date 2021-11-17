#include "../interface/Selection.h"

#include "../interface/Structures.h"

#define sqr(x) ((x)*(x))
#define min(x,y) (x<y ? x : y)

using namespace std;

const double nChi = 3.4; // FIXME from chi-distribution, 1 - exp(-x**2/2)

/*****************************************************************************/
Selection::Selection(int nCmTra) : nCmTra(nCmTra)
{
  const string pre = "../out/sel/"+to_string(nCmTra)+"part";

  // tracker
  his_cm_xy.init(-10,10,100, -10,10,100,  pre+"/cm_xy.his");
  his_cm_zr.init(-50,50,100,   0,20,100,  pre+"/cm_zr.his");

  his_cm_z.init(    -30,30,300,           pre+"/cm_z.his");
  his_cm_dz.init(   -30,30,300,           pre+"/cm_dz.his");
  his_cm_dzrel.init(-10,10,500,           pre+"/cm_dzrel.his");
  his_cm_dz_vs_sz.init(0,1,100, -2,2,100, pre+"/cm_dz_vs_sz.his");

  for(int i = 0; i < nParts; i++)
    his_cm_looper[i].init(0,2, 40, pre+"/cm_looper_"+parts[i]+".his", true); // normalize

  // roman pots
  his_rp_x12_elas.init(-0.1,0.1,100, -0.1,0.1,100,pre+"/rp_x12_elas.his");
  his_rp_x12_sign.init(-0.1,0.1,100, -0.1,0.1,100,pre+"/rp_x12_sign.his");
  his_rp_x12_side.init(-0.1,0.1,100, -0.1,0.1,100,pre+"/rp_x12_side.his");

  his_rp_dxrel.init(-50,50,100,                   pre+"/rp_dxrel.his");
  his_rp_dx_vs_sx.init(-0.1,0.1,100, -0.1,0.1,100,pre+"/rp_dx_vs_sx.his");

  // chi_2 vs chi_4
  his_chi_2.init( 0,20, 2000,          pre+"/chi_2.his" );
  his_chi_4.init( 0,20, 2000,          pre+"/chi_4.his" );
  his_chi_24.init(0,50, 250, 0,50,250, pre+"/chi_24.his");

  for(int topo = 0; topo < nTopos; topo++)
    his_chi_4_topo[topo].init( 0,20, 200,    pre+"/chi_4_"+topos[topo]+".his" );

  his_oth_2.init( 0,20, 200,           pre+"/oth_2.his" );
  his_oth_4.init( 0,20, 200,           pre+"/oth_4.his" );

  for(int iphi = 0; iphi < nDphi; iphi++)
  {
    string s = to_string(iphi);

    his_chi[iphi].init(0,20, 200, pre+"/chi_4_"+s+".his");
    his_oth[iphi].init(0,20, 200, pre+"/oth_4_"+s+".his");

    his_chi_24_dphi[iphi].init(0,50, 100, 0,50,100, pre+"/chi_24_"+s+".his");
  }

  // sump
  for(int topo = 0; topo < nTopos; topo++)
  {
    const string pre = "../out/rp/"+to_string(nCmTra)+"part";
    const string post = topos[topo]+".his";

    his_sum2p[topo].init(-0.5,0.5,200, -0.5,0.5,200, pre+"/sum2p_"+post);
    his_sum4p[topo].init(-0.5,0.5,200, -0.5,0.5,200, pre+"/sum4p_"+post);

    his_sumpx[topo].init(-1,1,200, -1,1,200, pre+"/sumpx_"+post);
    his_sumpy[topo].init(-1,1,200, -1,1,200, pre+"/sumpy_"+post);

    his_sum2px_vs_pred[topo].init(0,0.20,200, -0.5,0.5,200,
                                  pre+"/sum2px_vs_pred_"+post);
    his_sum2px_pred[topo].init(0,0.20,200,
                                  pre+"/sum2px_pred_"+post);
    his_sum2py[topo].init(-0.2,0.2,2000,
                                  pre+"/sum2py_"+post);
  }
}

/*****************************************************************************/
Selection::~Selection()
{
}

/*****************************************************************************/
bool Selection::areCmTracksLooper(const Event & event, int type, double w)
{
  double sump = (event.cmTracks[0].p +
                 event.cmTracks[1].p).length();

  his_cm_looper[type].fillw({sump / mass[type]}, w);

  return (sump / mass[type] < 0.2);
}

/*****************************************************************************/
bool Selection::areCmTracksPrimaries(const Event & event, double w)
{
  bool ok = true;

  // distance in z -> no cut
  const double dz = event.cmTracks[0].dz.val -
                    event.cmTracks[1].dz.val;

  double sz = sqrt(sqr(event.cmTracks[0].dz.sig) +
                   sqr(event.cmTracks[1].dz.sig));

  // z position -> cut
  for(auto & track : event.cmTracks)
  {
    // should come from near IP
    const double mean = -0.39;
    const double sig  =  4.59;

    // combined sigma
    double s = sqrt(sqr(sig) + sqr(track.dz.sig));

    if(fabs(track.dz.val - mean) > 4*s) // outside 4 sigma
      ok = false;
  }


  // from a photon conversion? -> cut
  const double & x = event.cmVertex.x.val;
  const double & y = event.cmVertex.y.val;
  const double & z = event.cmVertex.z.val;

  const double r = sqrt(sqr(x) + sqr(y));

  //
  his_cm_xy.fillw({x,y}, w);
  his_cm_zr.fillw({z,r}, w);

  if(r < 1) // primary if r < 1 cm
  {
    his_cm_z.fillw({event.cmTracks[0].dz.val}, w);
    his_cm_z.fillw({event.cmTracks[1].dz.val}, w);

          his_cm_dz.fillw({dz},    w);
       his_cm_dzrel.fillw({dz/sz}, w);
    his_cm_dz_vs_sz.fillw({sz,dz}, w);
  }
  else
    ok = false;

  return ok;
}

/*****************************************************************************/
bool Selection::areRpProtonsPrimaries(const Event & event, double w)
//                                      int evCategory)
{
  // mm -> cm
  const double x1 = event.prTracks[0].pos.x * 0.1; 
  const double x2 = event.prTracks[1].pos.x * 0.1;

  const double v1 = event.prTracks[0].pos.cxx * sqr(0.1);
  const double v2 = event.prTracks[1].pos.cxx * sqr(0.1);

  double dx = x1-x2;
  double sx = sqrt(v1+v2);

  if(nCmTra == 0)
  {
       his_rp_dxrel.fillw({dx/sx},w);
    his_rp_dx_vs_sx.fillw({sx,dx},w);
  }

  return (fabs(dx) < 0.01); // closer than 100 um FIXME VERY

/*
  if(evCategory == elastic ) his_rp_x12_elas.fillw({x1,x2},w); // elastic
  if(evCategory == signal  ) his_rp_x12_sign.fillw({x1,x2},w); // signal
  if(evCategory == sideband) his_rp_x12_side.fillw({x1,x2},w); // sideband

  if(evCategory == elastic)
  {
       his_rp_dxrel.fillw({dx/sx},w);
    his_rp_dx_vs_sx.fillw({sx,dx},w);
  }

  return (fabs(x1-x2) < 0.01); // closer than 100 um FIXME VERY
*/
}

/*****************************************************************************/
double Selection::getChi(const Vector2 & v)
{
  // https://en.wikipedia.org/wiki/Multivariate_normal_distribution
  // (cxx cxy)
  // (cxy cyy)

  double det = v.cxx * v.cyy - sqr(v.cxy);

  // x^T V^-1 x
  double xTVix = (v.cyy*sqr(v.x) - 2*v.cxy*v.x*v.y + v.cxx*sqr(v.y)) / det;

  // https://en.wikipedia.org/wiki/Chi_distribution
  // P(x) = x*exp(-x**2/2) [k=2]
  // CDF  = 1 - exp(-x**2/2)
  return sqrt(xTVix);
}

/*****************************************************************************/
void Selection::sumMomenta(const Event & event, Vector2 & sum2pt,
                                                Vector2 & sum4pt)
{
  sum2pt = {0,0, 0,0,0}; // 0-track or 2p   : x y cxx cxy cyy
  sum4pt = {0,0, 0,0,0}; // 2-track or 2p2h : x y cxx cxy cyy

  // collect roman pot protons
  for(auto & track : event.prTracks)
  {
    sum2pt += track.pt; // add proton
    sum4pt += track.pt; // add proton
  }

  // collect central hadrons
  for(auto & track : event.cmTracks)
  {
    // rotate variance of pt to variance of (px,py)
    double phi = track.p.phi();
    double c = cos(phi);
    double s = sin(phi);

    double V = sqr(track.sigpt);

    // set cm track.pt here, local
    Vector2 track_pt;

    track_pt.x = track.p.x; // copy
    track_pt.y = track.p.y;

    track_pt.cxx = c*c*V; // set up
    track_pt.cxy = c*s*V;
    track_pt.cyy = s*s*V;

    sum4pt += track_pt; // add hadron
  }
}

/*****************************************************************************/
int Selection::classifyEvent(const Event & event, int topo, double w)
{
  //
  double dphi = event.prTracks[0].p.dphi(
                event.prTracks[1].p);
  int iphi = int(dphi/(maxDphi/nDphi));
  if(iphi >= nDphi) iphi = nDphi-1;

  //
  vector<double> nTop = // FIXME ../out/sel/k_vs_dphi.par
    { 5.340, 5.339, 5.322, 5.307, 5.294, 5.267, 5.225, 5.196, 5.164,
      5.135, 5.137, 5.146, 5.176, 5.196, 5.233, 5.256, 5.286, 5.297 };

  ///////////////////////////////////////////////////////////
  // sum momenta
  Vector2 sum2pt, sum4pt;
  sumMomenta(event, sum2pt,sum4pt);

// modifying
/*
if(topo == TB) sum4pt.x -= -12e-3; // FIXME VERY
if(topo == BT) sum4pt.x -=   0e-3;
if(topo == TT) sum4pt.x -= -14e-3;
if(topo == BB) sum4pt.x -=   8e-3;
*/

  // collect
  his_sum2p[topo].fillw({sum2pt.x, sum2pt.y}, w);

  his_sum4p[topo].fillw({sum4pt.x, sum4pt.y}, w);

  his_sumpx[topo].fillw({sum2pt.x, sum4pt.x}, w);
  his_sumpy[topo].fillw({sum2pt.y, sum4pt.y}, w);

  his_sum2px_vs_pred[topo].fillw(
                         {sqrt(sum2pt.cxx), sum2pt.x}, w);

  his_sum2px_pred[topo].fillw(
                         {sqrt(sum2pt.cxx)}, w);

  if(nCmTra == 0) // FIXME VERY
    his_sum2py[topo].fillw({sum2pt.y}, w);
  else
    his_sum2py[topo].fillw({sum4pt.y}, w);

  // beam divergence
  // x ->  90e-6 * sqrt(2) / 45 * 6500e+3 = 18 MeV
  // y -> 100e-6 * sqrt(2) / 90 * 6500e+3 = 10 MeV

  // add beam divergence contribution
  const double sigmap_BD_x = 29e-3 * M_SQRT2; // 41.0 MeV // 29 vs 18 (pred) MeV
  const double sigmap_BD_y = 13e-3 * M_SQRT2; // 18.5 MeV // 13 vs 10 (pred) MeV

  sum2pt.cxx += sqr(sigmap_BD_x); sum2pt.cyy += sqr(sigmap_BD_y);
  sum4pt.cxx += sqr(sigmap_BD_x); sum4pt.cyy += sqr(sigmap_BD_y);

  ///////////////////////////////////////////////////////////
  // chi and event categories

  double chi_2p = getChi(sum2pt); // 0-track, p()p
  double chi_4p = getChi(sum4pt); // 2-track, p(h+h-)p

  //
  int evCategory = nothing;
 
  if(chi_2p < chi_4p)
  { // elasic
    if(chi_2p < nChi) evCategory = elastic;
  }
  else
  { // exclusive (signal or sideband)
    if(chi_4p < nChi) evCategory = signal;

    // \int_0^nChi x*exp(-k*x) = \int_nChi^nTop x*exp(-k*x)
    if(chi_4p > nChi &&
       chi_4p < nTop[iphi])
                      evCategory = sideband;
  }

  if(chi_4p > nChi)   his_chi_2.fillw({chi_2p},w);  // for elastic
  if(chi_2p > nChi) { his_chi_4.fillw({chi_4p},w);  // for exclusive
           his_chi_4_topo[topo].fillw({chi_4p},w);
                  his_chi[iphi].fillw({chi_4p},w); }

  if(evCategory == signal)    his_oth_2.fillw({chi_2p},w);
  if(evCategory == elastic) { his_oth_4.fillw({chi_4p},w);
                          his_oth[iphi].fillw({chi_4p},w); }

  his_chi_24.fillw({chi_2p,chi_4p},w);
  his_chi_24_dphi[iphi].fillw({chi_2p,chi_4p},w);

  return evCategory;
}

