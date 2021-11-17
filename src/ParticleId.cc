#include "../interface/ParticleId.h"

#include "../interface/MostProbable.h"
#include "../interface/Helper.h"
#include "../interface/Random.h"

#include "../interface/Structures.h"

#include <iostream>
#include <algorithm>

#define sqr(x) ((x)*(x))

using namespace std;

/*****************************************************************************/
// flag = 1 : collect
// flag = 2 : read and use relSigma

ParticleId::ParticleId(int flag)
{
  mostProbable = new MostProbable();

  if(flag == 1) // collect
  {
    // log(p), log(eps)
    elossAl.init(-3,3,150, 0.0,3.5,175, "../out/track/elossAl.his");
    elossNo.init(-3,3,150, 0.0,3.5,175, "../out/track/elossNo.his");
    elossPi.init(-3,3,150, 0.0,3.5,175, "../out/track/elossPi.his");
    elossKa.init(-3,3,150, 0.0,3.5,175, "../out/track/elossKa.his");
    elossPr.init(-3,3,150, 0.0,3.5,175, "../out/track/elossPr.his");
    elossSb.init(-3,3,150, 0.0,3.5,175, "../out/track/elossSb.his");

    for(int i = 0; i < nParts; i++)
    {
      Histo & h = relSigma[i];
      h.init(-maxEta,maxEta,nEta, 0,maxPt,nPt, 0.0,0.2,40,
             "../out/track/relSigma_"+parts[i]+".his");
    }
  }

  if(flag == 2) // for PID efficiency, for cmEff calculation
  {
    cerr << Helper::col(2) << " reading relsigmas for dE/dx" << Helper::col();

    for(int i = 0; i < nParts; i++)
    {
      Histo & h = relSigma[i];

      h.init(-maxEta,maxEta,nEta, 0,maxPt,nPt, 0.0,0.2,40,
             "../out/track/relSigma_"+parts[i]+".his");

      cerr << ".";
      h.read(); h.normalize();
    }

    cerr << " [done]" << endl;
  }
}

/*****************************************************************************/
ParticleId::~ParticleId()
{
  delete mostProbable;
}

/*****************************************************************************/
void ParticleId::processEnergyLoss(const vector<CmTrack> & tracks,
                                   int evCategory, int type)
{
  for(auto & tra : tracks)
  {
    //
    double eta = tra.p.eta();
    double pt  = tra.p.trans();

    double relsig = tra.sigma / tra.eps;

    if(type != unknown)
      relSigma[type].fill({eta,pt,relsig});

    //
    double logp   = log(tra.p.length());
    double logde  = log(tra.eps);

    if(evCategory == signal)
    {
                           elossAl.fill({logp,logde}); // all
      if(type == unknown)  elossNo.fill({logp,logde});
      if(type == pion)     elossPi.fill({logp,logde});
      if(type == kaon)     elossKa.fill({logp,logde});
      if(type == prot)     elossPr.fill({logp,logde});
    }

    if(evCategory == sideband)
                           elossSb.fill({logp,logde}); // sideband
  }
}

/*****************************************************************************/
// return particle type (unknown, pion, kaon, proton) for data
int ParticleId::identifyPair(const vector<CmTrack> & tra)
{
  const int nMulti = tra.size();

  if(tra.size() != 2)
  {
    cerr << " ParticleId::identifyPair problem: "
         << nMulti << " tracks" << endl;
    exit(1);
  }

  int type = unknown; // default

  vector<double> prob(nParts,1);

  for(int i = 0; i < nParts; i++)
  {
    const double & m = mass[i];

    vector<double> eps(nMulti);

    for(int j = 0; j < nMulti; j++)
      eps[j] = mostProbable->value(tra[j].p.length() / m);

    for(int j = 0; j < nMulti; j++)
      prob[i] *= exp(-sqr((tra[j].eps - eps[j])/tra[j].sigma)/2) / tra[j].sigma;
  }

  const float Cpid = 10; // PARAMETER

  if(prob[pion] > prob[kaon]*Cpid && prob[pion] > prob[prot]*Cpid) type = pion;
  if(prob[kaon] > prob[pion]*Cpid && prob[kaon] > prob[prot]*Cpid) type = kaon;
  if(prob[prot] > prob[pion]*Cpid && prob[prot] > prob[kaon]*Cpid) type = prot;

  return type;
}

/*****************************************************************************/
// return type //// isIdentified for simulation
int ParticleId::identifyPair(const vector<Vector3> & ps, int type)
{
  // generate dE/dx (eps and sigma) for the two tracks (ps)
  vector<CmTrack> tracks(2);

  for(int j = 0; j < 2; j++)
  { 
    tracks[j].p = ps[j];

    double eta = ps[j].eta();
    double pt  = ps[j].trans();
    if(pt >= 2) pt = 2 - 1e-3; // push below 2 GeV

    tracks[j].eps = mostProbable->value(ps[j].length() / mass[type]);

    double r = relSigma[type].sample({eta,pt}, theRandom);

    tracks[j].sigma = tracks[j].eps * r;

    tracks[j].eps += theRandom->getGauss() * tracks[j].sigma;
  }

  return identifyPair(tracks);
}

/*****************************************************************************/
void ParticleId::propeller(int i)
{
  char c[4] = {'|', '/','-','\\'};

  cerr << c[i % 4] << "\b";
}

/*****************************************************************************/
void ParticleId::prepareDemo(int type)
{
  const int nAll = 1e+5;
  const double dp = 50e-3;

  cerr << " simulating";

  ofstream fileOut("../out/track/pidDemo_"+parts[type]+".out");

  for(double p3_ = dp/2; p3_ < 2; p3_ += dp)
  {
  cerr << ".";

  for(double p4_ = dp/2;  p4_ < 2; p4_ += dp)
  {
    //
    vector<int> nId(nParts,0);

    for(int i = 0; i < nAll; i++) // PARS
    {
    double p3len = p3_ + theRandom->getFlat(-dp/2,dp/2);
    double p4len = p4_ + theRandom->getFlat(-dp/2,dp/2);

    double eta3 = theRandom->getFlat(-2.5,2.5);
    double eta4 = theRandom->getFlat(-2.5,2.5);

    double p3t = p3len / cosh(eta3);
    double p4t = p4len / cosh(eta4);

    double pz3 = p3t * sinh(eta3);
    Vector3 p3 = {p3t,0,pz3};

    double pz4 = p4t * sinh(eta4);
    Vector3 p4 = {p4t,0,pz4};

    vector<Vector3> ps; ps.push_back(p3); ps.push_back(p4);

    int pid = identifyPair(ps,type);

    if(pid != unknown) nId[pid]++;
    }

    fileOut << " " << p3_ << " " << p4_;
    for(int j = 0; j < nParts; j++)
      fileOut << " " << nId[j] / float(nAll);
    fileOut << endl;
  }
  fileOut << endl;
  }

  fileOut.close();

  cerr << " [done]" << endl;
}

/*****************************************************************************/
double ParticleId::getEfficiency(double eta3, double p3t,
                                 double eta4, double p4t, int type)
{
  double pz3 = p3t * sinh(eta3);
  Vector3 p3 = {p3t,0,pz3};

  double pz4 = p4t * sinh(eta4);
  Vector3 p4 = {p4t,0,pz4};

  vector<Vector3> ps; ps.push_back(p3); ps.push_back(p4);

  int nId = 0;
  for(int i = 0; i < nPidMcEvents; i++)
  {
    int pid = identifyPair(ps,type);

    if(pid == type) nId++;
  }

  return nId / float(nPidMcEvents);

/* one shot
 int pid = identifyPair(ps,type);

 if(pid == type) return 1.;
            else return 0.;
*/
}

