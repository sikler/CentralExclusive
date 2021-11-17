#include "../interface/RpReco.h"

#include "../interface/Structures.h"
#include "../interface/RpPat.h"
#include "../interface/RpFit.h"
#include "../interface/RpEffic.h"
#include "../interface/RpPolygon.h"

#include <iomanip>

#define sqr(x) ((x)*(x))

using namespace std;

/*****************************************************************************/
RpReco::RpReco(int nCmTra) : nCmTra(nCmTra)
{
  map<RpDet, map<vector<double>,int>> dummy;
  theRpFit   = new RpFit(dummy,nCmTra, 2); // read shifts, pattern fits
  theRpEffic = new RpEffic(2,nCmTra);      // read group effic

  // initialize optics [default numbers]

  // near
  v_x[0] = -2.24791387053766; // [mm]

  L_x[0][0] = 0.125396407127792E3;
  L_x[1][0] = 0.125396407127792E3;

  //
  v_y[0] = 0.025781593410852;

  L_y[0] = 238.517247191010E3;

  // far
  v_x[1] = -1.92610996810677;

  L_x[0][1] = -3.00655323980445E3;
  L_x[1][1] = -3.00655323980445E3;

  //
  v_y[1] = -0.000000021508565;

  L_y[1] = 271.511335947517E3;

  // overwrite using our covariance measuremen

  // arm 1
  {
    L_x[0][0] = 2660; // 1n
    L_x[0][1] = -830; // 1f
  }

  // arm 2
  {
    L_x[1][0] =   630; // 2n
    L_x[1][1] = -2580; // 2f
  }

  // initialize corrections [arm][sta][rpt]
  mx[0][1][0]=-0.465; mx[0][0][0]=-0.210;
  mx[1][0][0]= 0.167; mx[1][1][0]=-0.450;
  mx[0][1][1]=-0.081; mx[0][0][1]=-0.112;
  mx[1][0][1]= 0.373; mx[1][1][1]=-0.574;
  
  my[0][1][0]=-0.689; my[0][0][0]=-1.479;
  my[1][0][0]=-0.916; my[1][1][0]= 0.044;
  my[0][1][1]= 0.009; my[0][0][1]= 0.842;
  my[1][0][1]= 1.312; my[1][1][1]= 0.316;

  // additional shifts based on 0part TB and BT
  mx[0][0][0] -= 10e-3 - 8e-3; // 1nT
  mx[0][0][1] -= 41e-3; // 1nB
  mx[0][1][0] -= -5e-3 - 8e-3; // 1fT
  mx[0][1][1] -= 34e-3; // 1fB
  mx[1][0][0] -= -8e-3; // 2nT
  mx[1][0][1] -= 33e-3 - 8e-3; // 2nB
  mx[1][1][0] -= -9e-3; // 2fT
  mx[1][1][1] -= 31e-3 - 8e-3; // 2fB

  // additional shifts based on 0part sum p_y
  double dy;

  cerr << " dy [mm] :";

  dy = - 0.64*7000/(L_y[0]+L_y[1]); my[0][0][0]+=dy; my[0][1][0]-=dy; // 1nT 1fT
  cerr << " 1T=" << dy;
  dy = - 1.04*7000/(L_y[0]+L_y[1]); my[0][0][1]+=dy; my[0][1][1]-=dy; // 1nB 1fB
  cerr << " 1B=" << dy;
  dy = - 1.60*7000/(L_y[0]+L_y[1]); my[1][0][0]+=dy; my[1][1][0]-=dy; // 2nT 2fT
  cerr << " 2T=" << dy;
  dy = - 0.64*7000/(L_y[0]+L_y[1]); my[1][0][1]+=dy; my[1][1][1]-=dy; // 2nB 2fB
  cerr << " 2B=" << dy;

  cerr << endl;

  // histograms
  const string pre = "../out/rp/"+to_string(nCmTra)+"part";

  for(int arm = 0; arm < 2; arm++)
  for(int topo = 0; topo < nTopos; topo++)
  {
    const string post = to_string(arm)+"_"+topos[topo]+".his";

    his_prot_loc[arm][topo].init(-1,1,100, -4,4,100, pre+"/prot/loc_"+post);
    his_prot_mom[arm][topo].init(-1,1,100, -1,1,100, pre+"/prot/mom_"+post);
  }

  for(int topo = 0; topo < nTopos; topo++)
  {
    const string post = topos[topo]+".his";

    his_weight[topo].init(0,10,400, pre+"/weight_"+post);

    his_vtx_x[topo].init(-1,1,100, -1,1,100, pre+"/vtx/x_"+post);
    his_vtx_y[topo].init(-5,5,100, -5,5,100, pre+"/vtx/y_"+post);

    his_prot_mom_x[topo].init(-1,1,200, -1,1,200, pre+"/prot/mom_x_"+post);
    his_prot_mom_y[topo].init(-1,1,200, -1,1,200, pre+"/prot/mom_y_"+post);
  }

  for(int arm = 0; arm < 2; arm++)
  for(int sta = 0; sta < 2; sta++)
  for(int rpt = 0; rpt < 2; rpt++)
  for(int topo = 0; topo < nTopos; topo++)
  {
    const string post =
      to_string(arm)+to_string(sta)+to_string(rpt)+"_"+topos[topo]+".his";

    his_hits_loc[arm][sta][rpt][topo].init(-2,2,100, -40,40,100,
                                                    pre+"/hits/loc_"  +post);
    his_hits_loc_x[arm][sta][rpt][topo].init( -2, 2,100, pre+"/hits/loc_x_"+post);
    his_hits_loc_y[arm][sta][rpt][topo].init(-40,40,100, pre+"/hits/loc_y_"+post);
  }

  for(int arm = 0; arm < 2; arm++)
  for(int topo = 0; topo < nTopos; topo++)
  {
    const string post = to_string(arm)+"_"+topos[topo]+".his";

    his_hits_nf_x[arm][topo].init( -2, 8,200,  -2, 8,200, pre+"/hits/nf_x_"+post);
    his_hits_nf_y[arm][topo].init(-40,40,200, -40,40,200, pre+"/hits/nf_y_"+post);
  }

  //
  for(int a = 0; a < 2; a++)
  for(int s = 0; s < 2; s++)
  for(int r = 0; r < 2; r++)
  for(int u = 0; u < 2; u++)
  {
    RpDet det;
    det.arm = a; det.sta = s; det.rpt = r; det.uv  = u;

    char name[256];
    sprintf(name,"../out/rp/%dpart/occup/319311_%d%d%d_%d.his",
                 nCmTra, a,s,r, u);

    his_occup[det].init(0,nStrips,nStrips, -maxSlope,maxSlope,binSlope, name);
  }
}

/*****************************************************************************/
RpReco::~RpReco()
{
  ofstream fileTex("../an/rp_optics_"+to_string(nCmTra)+"part.tex");

  for(int arm = 0; arm < 2; arm++)
  for(int topo = 0; topo < nTopos; topo++)
  {
    double sig_xstar,sig_tstar, Ln,Lf;

    covariance(his_hits_nf_x[arm][topo], sig_xstar,sig_tstar, Ln,Lf);

    fileTex << fixed
      << " " << arm+1
      << " & " << topos[topo]
      << " & " << setprecision(0) << round(sig_xstar*1e+3)
      << " & " << setprecision(0) << round(sig_tstar*1e+6)
      << " & " << setprecision(0) << round(Ln/1e+1)*1e+1
      << " & " << setprecision(0) << round(Lf/1e+1)*1e+1
      << " & " << setprecision(arm==0?2:1) 
                                  << round(Lf/Ln*1e+2)/1e+2
      << " \\\\ " << endl;
  }

  fileTex.close();
}

/*****************************************************************************/
void RpReco::covariance(Histo & his,
                        double & sig_xstar, double & sig_tstar,
                        double & Ln, double & Lf)
{
  const double & v1 = v_x[0];
  const double & v2 = v_x[1];

  double c11=0, c12=0, c22=0, vt=0, sum=0;
  const double det = 7000.; // mm

  for(int ix = 0; ix < his.axes[0].bins; ix++)
  for(int iy = 0; iy < his.axes[1].bins; iy++)
  {
    double x = his.axes[0].center(ix);
    double y = his.axes[1].center(iy);

    if(abs(x) < 0.7 && abs(y) < 0.7) // [mm] // FIXME
    {
      double val = his.get({ix,iy}).first;

      c11 += x*x * val;
      c12 += x*y * val;
      c22 += y*y * val;

      vt += sqr((v2*x - v1*y)/det) * val; // var(theta*)

      sum += val;
    }
  }

  c11 /= sum;
  c12 /= sum;
  c22 /= sum;
  vt  /= sum;

  // var(x*)
  double vx = (c11*c22 - sqr(c12))/(c22*sqr(v1) - 2*c12*v1*v2 + c11*sqr(v2));

  sig_xstar = sqrt(vx);
  sig_tstar = sqrt(vt);

  double L1_ = sqrt( (c11 - vx*sqr(v1)) / vt );
  double L2_ = sqrt( (c22 - vx*sqr(v2)) / vt );

  Lf = -L2_;
  Ln =  L1_;
}

/*****************************************************************************/
bool RpReco::getTopo(const vector<RpTrack> & rpTracks, int & topo)
{
  int nConf = 0;
  topo = -1;

  // put together code
if(0)
  {
    char code[4] = {'x','x','x','x'};

    for(auto & track : rpTracks)
    {
      RpDet det = track.det;
   
      int loc = -1;

      if(det.arm == 0 && det.sta == 0) loc = 1;
      if(det.arm == 0 && det.sta == 1) loc = 0;

      if(det.arm == 1 && det.sta == 0) loc = 2;
      if(det.arm == 1 && det.sta == 1) loc = 3;

      code[loc] = (det.rpt == 0 ? 'T' : 'B');
    }

mtx_cout.lock();
    cout << " " << code << endl;
mtx_cout.unlock();
  }

  // must have all four tracklets
  if(rpTracks.size() == 4)
  {
    bool v[2][2][2];
    for(short int arm = 0; arm < 2; arm++)
    for(short int sta = 0; sta < 2; sta++)
    for(short int rpt = 0; rpt < 2; rpt++)
      v[arm][sta][rpt] = false;
    
    for(auto & track : rpTracks)
    {
      RpDet det = track.det;

      const short int & arm = det.arm;
      const short int & sta = det.sta;
      const short int & rpt = det.rpt;
      
      v[arm][sta][rpt] = true;
    }

    // trigger
    bool lT = (v[0][1][0] && v[0][0][0]);
    bool lB = (v[0][1][1] && v[0][0][1]);

    bool rT = (v[1][0][0] && v[1][1][0]);
    bool rB = (v[1][0][1] && v[1][1][1]);
    
    if(lT && rB) { topo = 0; nConf++; }
    if(lB && rT) { topo = 1; nConf++; }
    if(lT && rT) { topo = 2; nConf++; }
    if(lB && rB) { topo = 3; nConf++; }
  }
  
  return (nConf == 1);
}

/*****************************************************************************/
void RpReco::localReco(int run, int topo, RpTrack & track,
                       map<RpDet,double> & locCy,
                       map<RpDet,double> & locCx)
{
  RpDet det = track.det;

  vector<Vector1> lpos(2);
  vector<double> eff(2); // u/v

  for(int uv = 0; uv < 2; uv++)
  {
    det.uv = uv;

    // copy
    vector<double> strip;
    for(int pla = 0; pla < nPlanes; pla++)
      strip.push_back(track.clus[uv][pla]);

    vector<double> orig = strip;

    RpPat theRpPat(4, -1); // don't read, nCmTra is not used
    vector<int> base_step = theRpPat.toPattern(strip); // strip mod!

    //
    LocalFit localFit;

    if(theRpFit->hasPattern(det,strip))
    { // already fitted
      const int & ibas = base_step[0];
      const int & base = base_step[1];
      const int & step = base_step[2];

      localFit = theRpFit->getLocalPosition(det,strip, ibas,base,step);
    }
    else
    { // something new (why?), fit now
      mtx_rpfit.lock();
      localFit = theRpFit->fitTracklet(det,orig);
      mtx_rpfit.unlock();
    }

    // local hit position
    lpos[uv] = {localFit.Cy, sqrt(localFit.Vy)}; 

    // get efficiency
    eff[uv] = theRpEffic->getGroupEffic(run,det, localFit.Cy, localFit.Cx);

if(eff[uv] == 0)
{
  cerr << " " << run << " " << det.print()
       << " " << localFit.Cy
       << " " << localFit.Cx
       << endl;
  exit(1);
}

    locCy[det] = localFit.Cy;
    locCx[det] = localFit.Cx;
  } 

  track.weight = 1 / (eff[0] * eff[1]); // 1 / (tracklet uv[0] * tracklet uv[1])

  // get global hit position, overwrite!
  track.pos = theRpFit->getGlobalPosition(track.det, lpos);
}

/*****************************************************************************/
void RpReco::correctPos(RpTrack & track)
{
  const auto & det = track.det;

  const auto & arm = det.arm;
  const auto & sta = det.sta;
  const auto & rpt = det.rpt;

  track.pos.x += mx[arm][sta][rpt];
  track.pos.y += my[arm][sta][rpt];
}

/*****************************************************************************/
void RpReco::getStar(int arm, const Vector2 pos[], RpRes & res)
{
  // 0 = N, 1 = F

  // determinant, 7000 mm
  double Dx = v_x[0] * L_x[arm][1] -
              v_x[1] * L_x[arm][0];

  double Dy = v_y[0] * L_y[1] -
              v_y[1] * L_y[0];

  // theta*
  {
  const double a1 =   v_x[0]/Dx; // -0.000321131 
  const double a0 = - v_x[1]/Dx; //  0.000275159

  res.theta.x = a0*pos[0].x + a1*pos[1].x;

  const double b0 = 1/L_y[0]/2;  // 2.09628e-06
  const double b1 = 1/L_y[1]/2;  // 1.84154e-06

/*
  const vector<double> var = { pos[0].cyy, pos[1].cyy }; // FIXME

  const double av = 1/(1/var[0] + 1/var[1]);

  const double b0 = 1/L_y[0] * av/var[0];
  const double b1 = 1/L_y[1] * av/var[1];
*/

  res.theta.y = b0*pos[0].y + b1*pos[1].y;

  res.theta.cxx = sqr(a0)*pos[0].cxx + sqr(a1)*pos[1].cxx;
  res.theta.cyy = sqr(b0)*pos[0].cyy + sqr(b1)*pos[1].cyy;
  res.theta.cxy = a0* b0* pos[0].cxy + a1* b1* pos[1].cxy;
  }

  // x*, y*
  {
  const double a0 =   L_x[arm][1]/Dx; //  0.239 -0.566
  const double a1 = - L_x[arm][0]/Dx; // -0.798  0.141

  res.pos.x = a0*pos[0].x + a1*pos[1].x;

  const double b0 =   L_y[1]/Dy; //  34.0738924558586
  const double b1 = - L_y[0]/Dy; // -38.7873337067881

  res.pos.y = b0*pos[0].y + b1*pos[1].y;

  res.pos.cxx = sqr(a0)*pos[0].cxx + sqr(a1)*pos[1].cxx;
  res.pos.cyy = sqr(b0)*pos[0].cyy + sqr(b1)*pos[1].cyy;
  res.pos.cxy = a0* b0* pos[0].cxy + a1* b1* pos[1].cxy;
  }
}

/*****************************************************************************/
void RpReco::globalReco(int topo,
                        const vector<RpTrack> & rpTracks,
                              vector<PrTrack> & prTracks,
                        double & rpWeight)
{
  //
  Vector2 pos[2][2]; // [arm][sta]
  rpWeight = 1.;

  // each tracklet
  for(auto & track : rpTracks)
  {
    RpDet det = track.det;

    const short int & arm = det.arm;
    const short int & sta = det.sta;

    pos[arm][sta] = track.pos;

    // Roman pots weight = Prod(proton weights) = Prod(tracklet weights)
    rpWeight *= track.weight;
  }

  //
  for(auto & track : rpTracks)
  {
    const RpDet & det = track.det;

    const short int & arm = det.arm;
    const short int & sta = det.sta;
    const short int & rpt = det.rpt;

    his_hits_loc[arm][sta][rpt][topo].fillw({track.pos.x,
                                             track.pos.y}, rpWeight);

    his_hits_loc_x[arm][sta][rpt][topo].fillw({track.pos.x}, rpWeight);
    his_hits_loc_y[arm][sta][rpt][topo].fillw({track.pos.y}, rpWeight);
  }

  // nf
  for(int arm = 0; arm < 2; arm++)
  {
    his_hits_nf_x[arm][topo].fillw({pos[arm][0].x,
                                    pos[arm][1].x}, rpWeight);

    his_hits_nf_y[arm][topo].fillw({pos[arm][0].y,
                                    pos[arm][1].y}, rpWeight);
  }

  // fill weight
  his_weight[topo].fill({rpWeight});

  // each arm
  for(int arm = 0; arm < 2; arm++)
  {
    RpRes res;
    getStar(arm, pos[arm], res); // calculate

    // proton tracks
    PrTrack track;

       track.p.x = - BeamP      * res.theta.x; // p*
       track.p.y =   BeamP      * res.theta.y;

    // set rp track.pt here
    track.pt.x   = track.p.x; // copy
    track.pt.y   = track.p.y; // copy

    track.pt.cxx =   sqr(BeamP) * res.theta.cxx;
    track.pt.cxy = - sqr(BeamP) * res.theta.cxy;
    track.pt.cyy =   sqr(BeamP) * res.theta.cyy;

       track.p.z = (arm == 0 ? +BeamP : -BeamP);

       track.pos = res.pos; // x*, y*

    his_prot_loc[arm][topo].fillw({track.pos.x, track.pos.y}, rpWeight); 
    his_prot_mom[arm][topo].fillw({track.p.x  , track.p.y  }, rpWeight); 

    prTracks.push_back(track);
  }

  his_prot_mom_x[topo].fillw({prTracks[0].p.x, prTracks[1].p.x}, rpWeight); 
  his_prot_mom_y[topo].fillw({prTracks[0].p.y, prTracks[1].p.y}, rpWeight); 

  his_vtx_x[topo].fillw({prTracks[0].pos.x, prTracks[1].pos.x}, rpWeight); 
  his_vtx_y[topo].fillw({prTracks[0].pos.y, prTracks[1].pos.y}, rpWeight); 
}

/*****************************************************************************/
// reco and setp rpWeight
bool RpReco::reconstruct(int topo, Event & event)
{
  const int & run = event.run;
  vector<RpTrack> & rpTracks = event.rpTracks;
  vector<PrTrack> & prTracks = event.prTracks;

  // for his_occup
  map<RpDet,double> locCy,locCx;

  // re-reco and update rpTrack.pos
  for(auto & rpTrack : rpTracks)
  {
    localReco(run,topo, rpTrack, locCy,locCx);

    correctPos(rpTrack);
  }

  // check within local position is within limits
  bool isWithin = true;

  {
    const double scale = 271.5/238.5; // ~ 1.14 near -> far

    for(auto & rpTrack : rpTracks)
    {
      const double y_low =  7 * (rpTrack.det.sta == 0 ? 1 : scale); // n f
      const double y_hig = 24 * (rpTrack.det.sta == 0 ? 1 : scale); // [mm]

      if(fabs(rpTrack.pos.y) < y_low ||
         fabs(rpTrack.pos.y) > y_hig) isWithin = false;
    }
  }

  if(isWithin)
  {
    // global reco
    mtx_his.lock();
    globalReco(topo, rpTracks, prTracks, event.rpWeight); // get rpWeight!
    mtx_his.unlock();
  }

  // for his_occup
  if(isWithin)
  if(topo == TT || topo == BB)
  for(auto & rpTrack : rpTracks)
  {
    RpDet det = rpTrack.det;

    for(int uv = 0; uv < 2; uv++)
    {
      det.uv = uv;

      mtx_his.lock();
      his_occup[det].fillw({locCy[det], locCx[det]}, event.rpWeight);
      mtx_his.unlock();
    }
  }

  return isWithin;
}

