#include "../interface/Physics.h"

#include "../interface/Structures.h"

#include "../interface/RpVeto.h"
#include "../interface/TrkReco.h"
#include "../interface/Distributions.h"
#include "../interface/Resolutions.h"

#include "../interface/Helper.h"

#include <iostream>

using namespace std;

#define sqr(x) ((x)*(x))

/*****************************************************************************/
Physics::Physics(int nCmTra)
{
  theRpVeto  = new RpVeto(2, nCmTra); // read angular coverage
  theTrkReco = new TrkReco(2);        // read and use cmEff, and sim_effOne

  theDistributions = new Distributions();
  theResolutions   = new Resolutions(2); // read pt reso, collect mass resp

  // read effective integrated lumi per topo [(\int L) exp(-mu)]
  cerr << Helper::col(2) << " reading effective integrated lumi.."
       << Helper::col();

  string s;
  ifstream file("../out/lumi/int_lumi.dat");

  for(int topo = 0; topo < nTopos; topo++)
  {
    file >> s                      // topo
         >> s >> s                 // & recLumi
         >> s >> effIntLumi[topo]; // & effLumi

    // convert b^{-1} -> μb^{-1}
    effIntLumi[topo] *= 1e-6;
  }
  file.close();

  cerr << " [done]" << endl;
}

/*****************************************************************************/
Physics::~Physics()
{
  delete theDistributions;
  delete theResolutions;
}

/*****************************************************************************/
void Physics::correctMomenta(const Event & event,
                             Vector4 & pa, Vector4 & pb,
                             Vector4 & p1, Vector4 & p2, vector<Vector4> & hs)
{
  pa = Vector4::fourVector(mass[prot], {0,0,+BeamP});
  pb = Vector4::fourVector(mass[prot], {0,0,-BeamP});

  p1 = Vector4::fourVector(mass[prot], event.prTracks[0].p);
  p2 = Vector4::fourVector(mass[prot], event.prTracks[1].p);

  double E = 0, pz = 0;
  for(auto & track : event.cmTracks)
  {
    hs.push_back(Vector4::fourVector(mass[event.type], track.p));

    E  += hs.back().E;
    pz += hs.back().p.z;
  }

  // use energy-momentum conservation to correct longi momenta of protons
  p1.p.z += - E/2 - pz/2;
  p2.p.z +=   E/2 - pz/2;

  // do not use momentum conservation to correct (transverse) momenta of protons

  // recalculate energies for protons
  p1.E = sqrt(p1.p.length2() + sqr(mass[prot]));
  p2.E = sqrt(p2.p.length2() + sqr(mass[prot]));
}

/*****************************************************************************/
// get angles based on axes
void Physics::getAngles
  (const Vector3 & ay, const Vector3 & az,
   const Vector3 & p, vector<double> & a)
{
  const Vector3 ax = ay ^ az;

  double cosTheta =       az * p.norm();
  double phi      = atan2(ay * p, ax * p);

  a = {cosTheta,phi};
}

/*****************************************************************************/
// calculate angles for GJ frames
void Physics::calculateAngles(const Vector4 & pa, const Vector4 & pb,
                                    Vector4 q1, Vector4 q2,
                                    Vector4 p3, Vector4 p4,
                                    vector<double> & a)
{
  Boost b_to_cm = Vector4::getBoost(p3 + p4);

  Vector3 ay,az;

  if(fabs(q1.mass2()) > fabs(q2.mass2()) )
  {
    ay = (q1.p ^ q2.p).norm();
    p3 ^ b_to_cm; q1 ^ b_to_cm; q2 ^ b_to_cm;
    az = q1.p.norm();
  }
  else
  {
    ay = (q2.p ^ q1.p).norm();
    p3 ^ b_to_cm; q2 ^ b_to_cm; q1 ^ b_to_cm;
    az = q2.p.norm();
  }

  //
  getAngles(ay,az, p3.p, a);
}

/*****************************************************************************/
double Physics::getPhi_tilde(Vector4 p1, Vector4 p2, Vector4 q1,
                             const Boost & b_to_cm)
{
  // boost
  p1 ^ b_to_cm; p2 ^ b_to_cm;
  q1 ^ b_to_cm; 
  
  Vector3 a = (p1.p ^ q1.p);
  Vector3 b = (p2.p ^ q1.p);
   
  const double dphi_tilde = acos( (a * b) / a.length() / b.length() );
     
  return dphi_tilde;
} 

/*****************************************************************************/
void Physics::processMomenta(const Event & event,
                             const Vector4 & pa, const Vector4 & pb,
                             const Vector4 & p1, const Vector4 & p2,
                             const vector<Vector4> & hs)
{   
  //
  const Vector4 & p3 = hs[0];
  const Vector4 & p4 = hs[1];
    
  // resonance
  Vector4 p = p3 + p4;
      
  double pt = p.p.trans();
  double m  = p.mass();

  //
  Vector4 q1 = pa - p1;
  Vector4 q2 = pb - p2;

  vector<double> a;
  calculateAngles(pa,pb, q1,q2, p3,p4, a);

  const double dphi = fabs(p1.p.dphi(p2.p));
  
  // dphi_tilde
  const Boost b_to_cm = Vector4::getBoost(p3+p4);
  const double dphi_tilde = getPhi_tilde(p1,p2, q1, b_to_cm);
  
  //
  vector<float> dphip1tp2t = {dphi, p1.p.trans(), p2.p.trans()};
    
  vector<float> ptmcosthetaphi = {pt,m, a[0],a[1]};
 
  // do physics
  if(fabs(p.rap()) < maxResoY) // PARAMETER FIXME
//  if(fabs(p3.rap()) < maxHadronY &&
//     fabs(p4.rap()) < maxHadronY)
  {
//  const float & dphi = dphip1tp2t[0];
    const float & p1t  = dphip1tp2t[1];
    const float & p2t  = dphip1tp2t[2];

    const float & m   = ptmcosthetaphi[1];

    if(p1t > minKt && p1t < maxKt &&
       p2t > minKt && p2t < maxKt &&
       m < maxMass)
    if(theTrkReco->notFake(p3.p, p4.p, event.type))
    {
      // roman pot track weight; calculated in RpReco, 1/(prod of tracklet eff)
      double  rpWeight = event.rpWeight;

      double  rpCover  = theRpVeto->getAngularCoverage(event.topo, dphip1tp2t);

      double  cmEffic  = theTrkReco->getCmEff(dphip1tp2t,ptmcosthetaphi,
                                              event.topo, event.type);

      // calculate weight 
      double weight = event.sign / effIntLumi[event.topo]
                    * (rpWeight / rpCover) / cmEffic; // FIXME

      // for dphi and mass
if(cmEffic > 1e-2) // FIXME VERY
{
      theDistributions->collectForPhysics(event.topo, event.type,
                                          dphi_tilde, dphip1tp2t, m,
                                          q1,q2, p3,p4, weight);

      // for mass resolutions
      if(event.cat == signal)
        theResolutions->collectForMass(event.type, dphip1tp2t, m, 
                                       rpWeight / cmEffic, p3,p4); // FIXME
}
      //
    }
  }
} 

/*****************************************************************************/
void Physics::process(const Event & event)
{
  // correct and get momenta
  Vector4 pa,pb,p1,p2;
  vector<Vector4> hs;
  correctMomenta(event, pa,pb,p1,p2,hs);

  // process momenta
  processMomenta(event, pa,pb,p1,p2,hs);
}
