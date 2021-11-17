#include "../interface/TrkReco.h"

#include "../interface/TrkEffic.h"
#include "../interface/ParticleId.h"
#include "../interface/Random.h"
#include "../interface/RandomSobol.h"

#include "../interface/Helper.h"
#include "../interface/gzstream.h"

#include "../interface/Structures.h"

#include <thread> // for reading efficiencies only

#define sqr(x) ((x)*(x))

using namespace std;

/*****************************************************************************/
// flag = 1 : collect
// flag = 0 : calclate cmEff
// flag = 2 : read and use cmEff

TrkReco::TrkReco(int flag)
{
  if(flag == 1) // collect [data]
  {
     his_trackZ.init(-50,50, 1000, "../out/track/trk_z.his", true); // z of tra
     his_trackT.init( -5, 5, 1000, "../out/track/trk_t.his", true); // dT of tra
    his_trackDz.init(-20,20, 1000, "../out/track/trk_dz.his",true); // z4 - z3

    his_trackEtaPtPos.init(-maxEta,maxEta,nEta,
                                 0,maxPt ,mPt ,"../out/track/etaPtPos.his");
    his_trackEtaPtNeg.init(-maxEta,maxEta,nEta,
                                 0,maxPt ,mPt ,"../out/track/etaPtNeg.his");

    his_trackEtaPhiPos.init(-maxEta,maxEta,nEta,
                            -maxPsi,maxPsi,nPsi,
                            "../out/track/etaPhiPos.his");

    his_trackEtaPhiNeg.init(-maxEta,maxEta,nEta,
                            -maxPsi,maxPsi,nPsi,
                            "../out/track/etaPhiNeg.his");

    //
    for(int i = 0; i < nParts; i++)
    {
/*
      { Histo & h = his_trackSp[i];
        h.init(0,2, 40, "../out/track/trk_sp_"+parts[i]+".his", true); }
*/

      { Histo & h = his_trackChi2[i];
        h.init(-0.5,49.5,50, 0,100,200,
                      "../out/track/trk_chi2_"+parts[i]+".his", true); }

      { Histo & h = his_trackY[i];
        h.init(-5,5,200,
                         "../out/track/trk_y_"+parts[i]+".his", true); }

      { Histo & h = his_resoY[i];
        h.init(-5,5,200,
                         "../out/track/reso_y_"+parts[i]+".his", true); }
    }
  }

  if(flag == 2) // physics [data]
  {
    readCmEffTables();
  }

  if(flag == 0 || flag == 2) // for cmEff and for physics
  {
    //
    cerr << Helper::col(2) << " reading one-track efficiencies"
         << Helper::col();

    for(int type = 0; type < nParts; type++)
    {
      cerr << ".";

      //
      his_oneTrack_pos[type].init(-maxEta,maxEta,nEta, 0,maxPt,mPt,
                            "../out/track/sim_effOne_"+parts[type]+"p.his");
      his_oneTrack_pos[type].read();

      //
      his_oneTrack_neg[type].init(-maxEta,maxEta,nEta, 0,maxPt,mPt,
                            "../out/track/sim_effOne_"+parts[type]+"m.his");
      his_oneTrack_neg[type].read();
    }

    cerr << " [done]" << endl;
  }
}

/*****************************************************************************/
TrkReco::~TrkReco()
{
}

/*****************************************************************************/
void TrkReco::collectPars(const vector<CmTrack> & tracks)
{
  const double maxDt   =  2; // cm
  const double maxDz   = 15; // cm
  const double maxDiff = 10; // cm

  double diff = tracks[1].dz.val - tracks[0].dz.val;

  for(size_t j = 0; j < tracks.size(); j++)
  {
    if(fabs(tracks[j].dt.val)<maxDt && fabs(diff)<maxDiff)
       his_trackZ.fill({tracks[j].dz.val});

    if(fabs(tracks[j].dz.val)<maxDz && fabs(diff)<maxDiff)
       his_trackT.fill({tracks[j].dt.val});
  }

  if(fabs(tracks[0].dt.val) < maxDt && fabs(tracks[0].dz.val) < maxDz &&
     fabs(tracks[1].dt.val) < maxDt && fabs(tracks[1].dz.val) < maxDz)
    his_trackDz.fill({diff});

  for(int i = 0; i < 2; i++)
  {
    double eta = tracks[i].p.eta();
    double  pt = tracks[i].p.trans();

    if(tracks[i].q == 1) his_trackEtaPtPos.fill({eta,pt});
                    else his_trackEtaPtNeg.fill({eta,pt});

    double phi = tracks[i].p.phi();

    if(tracks[i].q == 1) his_trackEtaPhiPos.fill({eta,phi});
                    else his_trackEtaPhiNeg.fill({eta,phi});
  }
}

/*****************************************************************************/
void TrkReco::collectChi2(const vector<CmTrack> & tracks, int type)
{
  for(size_t j = 0; j < tracks.size(); j++)
  {
    //
    his_trackChi2[type].fill({float(tracks[j].ndf), tracks[j].chi2});

    //
    Vector4 p = Vector4::fourVector(mass[type], tracks[j].p);

    his_trackY[type].fill( { p.rap() } );
  }


  Vector4 p3 = Vector4::fourVector(mass[type], tracks[0].p);
  Vector4 p4 = Vector4::fourVector(mass[type], tracks[1].p);
  double y = (p3+p4).rap();

  his_resoY[type].fill({y});
}

/*****************************************************************************/
bool TrkReco::hasGoodEffic(const Vector3 & p)
{
  return (fabs(p.eta())  < 2.50 &&  // |eta| < 2.5  // FIXME
               p.trans() > 0.10);   // pt > 0.1 GeV // FIXME
}

/*****************************************************************************/
bool TrkReco::hasGoodEffic(const vector<CmTrack> & tracks)
{
  return (hasGoodEffic(tracks[0].p) &&
          hasGoodEffic(tracks[1].p)); 
}

/*****************************************************************************/
// input:
// incoming
//  pa = Vector4::fourVector(mass[prot], {0,0,+BeamP});
//  pb = Vector4::fourVector(mass[prot], {0,0,-BeamP});
// scattered protons
//  p1 = Vector4::fourVector(mass[prot], {p1x,p1y, 0});
//  p2 = Vector4::fourVector(mass[prot], {p2x,p2y, 0});

void TrkReco::getMomenta(const Vector4 & pa, const Vector4 & pb,
                               Vector4   p1,       Vector4   p2,
                          int type, const vector<double> & ptmcosthetaphi,
                          double y,
                               Vector4 & p3,       Vector4 & p4)
{
  const double & pt  =      ptmcosthetaphi[0];
  const double & m   =      ptmcosthetaphi[1];
  const double theta = acos(ptmcosthetaphi[2]);
  const double & phi =      ptmcosthetaphi[3];

  const double mt = sqrt(sqr(m) + sqr(pt));

  const double deltap1z = -mt/2 * exp( y);
  const double deltap2z =  mt/2 * exp(-y);

  // set pz for scattered protons
  p1.p.z = pa.p.z + deltap1z;
  p2.p.z = pb.p.z + deltap2z;

  //
        Vector4 q1 = pa - p1;
  const Vector4 q2 = pb - p2;

  // GJ
  Vector3 ay = (q1.p ^ q2.p).norm();

  // boost to cm
  const Vector3 p = (q1 + q2).p;
  const double  E = sqrt(sqr(m) + p.length2());

  const double beta = p.length() / E;
  const Vector3 n   = p.norm();

  const Boost b_to_cm = pair<double,Vector3>(beta,n);

  //
  q1 ^ b_to_cm;

  //
  const Vector3 az = q1.p.norm();
  const Vector3 ax = ay ^ az;

  //
  double q = sqrt(sqr(m/2) - sqr(mass[type]));

  const Vector3 dir = (ax * sin(theta)*cos(phi)) +
                      (ay * sin(theta)*sin(phi)) +
                      (az * cos(theta));

  const Vector3 mom = dir * q;
  
  p3 = Vector4::fourVector(mass[type], mom * (+1));
  p4 = Vector4::fourVector(mass[type], mom * (-1));
    
  // boost p3,p4 to lab
  const Boost b_to_lab = pair<double,Vector3>(-beta,n);
  p3 ^ b_to_lab; p4 ^ b_to_lab;
}

/*****************************************************************************/
bool TrkReco::notFake(const Vector3 & p3,
                      const Vector3 & p4, int type)
{
  double oneTrackEff_pos = his_oneTrack_pos[type].val({p3.eta(),p3.trans()});
  double oneTrackEff_neg = his_oneTrack_neg[type].val({p4.eta(),p4.trans()});

  return(oneTrackEff_pos > minTkEff &&
         oneTrackEff_neg > minTkEff);
}

/*****************************************************************************/
// -cmEff
void TrkReco::prepareCmTables(int type)
{
  Random theRandom;

  TrkEffic theTrkEffic(1);     // read and use sim tables
  ParticleId theParticleId(2); // read and use relSigma

  const int longSobol = 256; // per bin FIXME
  const int nSobol = 8;      // number of elements in random vector
  RandomSobol theRandomSobol(longSobol, nSobol);

  //
  cerr << Helper::col(1) << " calculating combined tracker"
       << " [" << nDphi<<"*"<<nMass<<"*"<<nCosTheta<<"*"<<nPhi << " bins]"
       << " [" << longSobol << " per bin]"
       << " * [" << nPidMcEvents << " for Pid]"
       << Helper::col() << endl;

  const double dKt = (maxKt - minKt)/nKt;

  //
  for(int i1 = 0;  i1 < nKt; i1++)
  for(int i2 = i1; i2 < nKt; i2++)
  {
    HistoVal his_cm_eff[nTopos], his_cm_all[nTopos];
//    Histo his_cm_eff[nTopos], his_cm_all[nTopos];

    double low1 = minKt + i1*dKt;
    double low2 = minKt + i2*dKt;

    double hig1 = minKt + (i1+1)*dKt;
    double hig2 = minKt + (i2+1)*dKt;

    for(int topo = 0; topo < nTopos; topo++)
    {
      his_cm_eff[topo].init(0,maxDphi,nDphi, // dphi
                            0,maxMass,nMass, // m
                            -maxCosTheta, maxCosTheta, nCosTheta, // cos(theta)
                            -maxPhi,      maxPhi,      nPhi,      // phi
   "../out/track/cmEff/"+parts[type]+"_"+to_string(i1)
                                    +"_"+to_string(i2)+"_"+topos[topo]+".his");

      his_cm_all[topo].init(0,maxDphi,nDphi, // dphi
                            0,maxMass,nMass, // m
                            -maxCosTheta, maxCosTheta, nCosTheta, // cos(theta)
                            -maxPhi,      maxPhi,      nPhi,      // phi
                            "");
    }

    cerr << " (" << i1 << "," << i2 << ") ";

    const double _dphi = maxDphi/nDphi; // [0,pi]
    const double _mass = maxMass/nMass; // [0,4]
    const double _cthe = 2*maxCosTheta/nCosTheta; // [-1,1]
    const double _phi  = 2*maxPhi     /nPhi;      // [-pi:pi]

    const int imass0 = int(2*mass[type] / _mass);

    int iprop = 0, jprop = 0;

    //
    for(int idphi = 0; idphi < nDphi; idphi++)
    for(int imass = imass0;
                       imass < nMass; imass++)
    for(int icthe = 0; icthe < nCosTheta; icthe++)
    for(int iphi  = 0; iphi  < nPhi;      iphi++)
    {
      if(++iprop % (nCosTheta*nPhi) == 0)
      {
        if(iprop % (nMass*nCosTheta*nPhi) == 0) cerr << ".";
        else
        {
          char c[4] = {'|', '/','-','\\'};
          cerr << c[jprop] << "\b";

          jprop = (jprop+1) % 4;
        }
      }

      //
      theRandomSobol.reset();

      for(int i = 0; i < longSobol; i++)
      {
      vector<double> r = theRandomSobol.get();
      int ir = 0;

      // scattered protons in roman pots, using 0,1
      double p1t  = theRandom.getFlat(low1,hig1, r[ir++]);
      double p2t  = theRandom.getFlat(low2,hig2, r[ir++]);

      // using 2,3,4,5
      double dphi = idphi*_dphi + theRandom.getFlat(0,_dphi, r[ir++]);
      double m    = imass*_mass + theRandom.getFlat(0,_mass, r[ir++]);

      if(m < 2*mass[type]) continue;

      double cosTheta =
                    -maxCosTheta +
                              icthe*_cthe + theRandom.getFlat(0,_cthe, r[ir++]);
      double phi  = -maxPhi + iphi *_phi  + theRandom.getFlat(0,_phi , r[ir++]);

      double phi1 = theRandom.getFlat(0,2*M_PI,  r[ir++]); // using 6
      double phi2 = phi1 + dphi;

      //
      double p1x = p1t * cos(phi1);
      double p1y = p1t * sin(phi1);

      double p2x = p2t * cos(phi2);
      double p2y = p2t * sin(phi2);

      // get topology
      int topo = -1;
      if(p1y >= 0 && p2y <  0) topo = TB;
      if(p1y <  0 && p2y >= 0) topo = BT;
      if(p1y >= 0 && p2y >= 0) topo = TT;
      if(p1y <  0 && p2y <  0) topo = BB;
      if(topo == -1) exit(1);

      // resonance
      double pt = sqrt(sqr(p1x + p2x) + sqr(p1y + p2y));

      // max rapidity FIXME here reso rap!
      double y = maxResoY * theRandom.getFlat(-1,1, r[ir++]); // using 7

      // incoming protons
      Vector4 pa,pb, p1,p2;
      pa = Vector4::fourVector(mass[prot], {0,0,+BeamP});
      pb = Vector4::fourVector(mass[prot], {0,0,-BeamP});

      // scattered protons
      p1 = Vector4::fourVector(mass[prot], {p1x,p1y, 0});
      p2 = Vector4::fourVector(mass[prot], {p2x,p2y, 0});

      Vector4 p3,p4; // pos,neg
      getMomenta(pa,pb, p1,p2, type,{pt,m,cosTheta,phi},y, p3,p4);

      //
      double trkEff = 0; 
      // get from database
      if(notFake(p3.p, p4.p, type)) // FIXME
        trkEff = theTrkEffic.getEfficiency(p3.p, p4.p, type);

      // MC calculation with nPidMcEvents=1 trials
      double pidEff =
              theParticleId.getEfficiency(p3.p.eta(), p3.p.trans(),
                                          p4.p.eta(), p4.p.trans(), type);
      //
      double w = trkEff * pidEff;

      his_cm_eff[topo].fillw({dphi,m, cosTheta,phi}, w);
      his_cm_all[topo].fill( {dphi,m, cosTheta,phi});
      }
    }

    for(int topo = 0; topo < nTopos; topo++)
      his_cm_eff[topo].div(his_cm_all[topo]);

    cerr << " [done]" << endl;
  }
}

/*****************************************************************************/
void TrkReco::readCmEff(const int type)
{
  mtx_cerr.lock();
  cerr << "(" << parts[type] << ")";
  mtx_cerr.unlock();

  for(int topo = 0; topo < nTopos; topo++)
  {
  for(int i1 = 0;  i1 < nKt; i1++)
  for(int i2 = i1; i2 < nKt; i2++)
  {
    HistoVal & h = his_cm_eff[topo][type][i1][i2];
//    Histo & h = his_cm_eff[topo][type][i1][i2];

    h.init(0,maxDphi,nDphi, // dphi
           0,maxMass,nMass, // m
           -maxCosTheta, maxCosTheta, nCosTheta, // cos(theta)
           -maxPhi,      maxPhi,      nPhi,      // phi
           "../out/track/cmEff/"+parts[type]+"_"+to_string(i1)
                                            +"_"+to_string(i2)
           +"_"+topos[topo]+".his.gz");

    h.read(true); // fromZip
  }

  mtx_cerr.lock();
  cerr << ".";
  mtx_cerr.unlock();
  }
}

/*****************************************************************************/
void TrkReco::readCmEffTables()
{
  int list;

  cerr << Helper::col(2) << " reading combined tracker efficiency"
       << Helper::col();

  list = pion; thread thpi(&TrkReco::readCmEff,this,list);
  list = kaon; thread thka(&TrkReco::readCmEff,this,list);
  list = prot; thread thpr(&TrkReco::readCmEff,this,list);

  thpi.join();
  thka.join();
  thpr.join();

  cerr << " [done]";
}

/*****************************************************************************/
double TrkReco::getCmEff(const vector<float> & dphip1tp2t,
                         const vector<float> & ptmcosthetaphi,
                         int topo, int type)
{
  const float & dphi = dphip1tp2t[0];
        float   p1t  = dphip1tp2t[1];
        float   p2t  = dphip1tp2t[2];

  // exchange, force p1t < p2t
  if(p1t > p2t)
  { double v = p1t; p1t = p2t; p2t = v; }

  const double & m        = ptmcosthetaphi[1];
  const double & costheta = ptmcosthetaphi[2];
  const double & phi      = ptmcosthetaphi[3];

  const double dKt = (maxKt - minKt)/nKt;

  //
  const int i1 = int((p1t - minKt)/dKt);
  const int i2 = int((p2t - minKt)/dKt);

  return his_cm_eff[topo][type][i1][i2].val({dphi,m, costheta,phi});
}

