#include "../interface/TrkEffic.h"

#include "../interface/Structures.h"

#include "../interface/Resolutions.h"
#include "../interface/Random.h"

#include "../interface/Helper.h"
#include "../interface/gzstream.h"

#include <thread> // for reading efficiencies only

#define sqr(x) ((x)*(x))

using namespace std;

//
//
const int minLayers   = 3;
const int minClusters = 5;

const int minBpixTracks  = 1;
const int minPixelTracks = 1;

const int mPt = (5 * nPt)/2; // FIXME

/*****************************************************************************/
// flag = 0 : process simulation
// flag = 1 : read and use sim tables

TrkEffic::TrkEffic(int flag)
{
  if(flag == 1)
    readSimTables();
}

/*****************************************************************************/
TrkEffic::~TrkEffic()
{
}

/*****************************************************************************/
struct Part
{
  int pid;
  float eta,pt,phi, vz;
  int nMatched;

  void print(const string & s)
  {
    cerr << " " << s << " " << pid
         << " (" << eta << " " << pt << " " << phi << ")"
         << " " << nMatched << endl;
  }
};

/*****************************************************************************/
int TrkEffic::getIndex(double low, double hig, int nbins, double x)
{
  int ix;

  if(x >= low && x < hig)
  {
    ix = int((x - low)/(hig - low) * nbins);

    if(ix < 0)      ix = 0;
    if(ix >= nbins) ix = nbins - 1;
  }
  else
    ix = -1;

  return ix;
}

/*****************************************************************************/
bool TrkEffic::getIndex(const Part & p, vector<short int> & index)
{
  short int ieta = getIndex(-maxEta,maxEta, nEta, p.eta);
  short int ipt  = getIndex(      0,maxPt , mPt , p.pt );
  short int ipsi = getIndex(-maxPsi,maxPsi, nPsi, p.phi);

  index = {ieta,ipt,ipsi};

  return (!(ieta==-1 || ipt==-1 || ipsi==-1));
}

//
struct SimTrackInfo {
  int nall;
  int nrec;
  int nhlt;
  map<pair<unsigned char,unsigned char>,int> noHltCodes;
};

/*****************************************************************************/
int TrkEffic::getLayerCountBpix(int layerCodes)
{
  int count = 0;

  for(int i = 0; i < 4; i++)
    count += (layerCodes & (1<<i) ? 1 : 0);

  return count;
}

/*****************************************************************************/
bool TrkEffic::isCompatible(const Part & gen, const Part & rec)
{
  if(rec.pid * gen.pid > 0) // same sign
  {
    vector<double> g = {gen.pt * cos( gen.phi),
                        gen.pt * sin( gen.phi),
                        gen.pt * sinh(gen.eta)};

    vector<double> r = {rec.pt * cos( rec.phi),
                        rec.pt * sin( rec.phi),
                        rec.pt * sinh(rec.eta)};

    double dif = sqrt(sqr(g[0] - r[0]) +
                      sqr(g[1] - r[1]) +
                      sqr(g[2] - r[2]));

    double len = gen.pt * cosh(gen.eta);

    if(dif < 0.25 * len) // FIXME absolute resolution better than 25%
      return true;
  }

  return false;
}

/*****************************************************************************/
void TrkEffic::processSimulation(const string & partName)
{
  cerr << Helper::col(1) << " processing simulation" << Helper::col() << endl;

  map<vector<short int>,SimTrackInfo> simInfo; // local

  //
  int nEvents = 0;

  Histo his_genTrack, his_oneTrack, his_mulTrack, his_hltTrack;

  Histo his_genEtaPhi, his_oneEtaPhi;
  Histo his_genPtPhi, his_onePtPhi;

  his_genTrack.init(-maxEta,maxEta,nEta, 0,maxPt,mPt, "");
  his_oneTrack.init(-maxEta,maxEta,nEta, 0,maxPt,mPt,
                  "../out/track/sim_effOne_"+partName+".his");
  his_mulTrack.init(-maxEta,maxEta,nEta, 0,maxPt,mPt,
                  "../out/track/sim_effMul_"+partName+".his");

  his_hltTrack.init(-maxEta,maxEta,nEta, 0,maxPt,mPt,
                  "../out/track/sim_effHlt_"+partName+".his");

  his_genEtaPhi.init(-maxEta,maxEta,nEta, -maxPsi,maxPsi,nPsi, "");
  his_oneEtaPhi.init(-maxEta,maxEta,nEta, -maxPsi,maxPsi,nPsi,
                  "../out/track/sim_effOne_etaPhi_"+partName+".his"); 

  his_genPtPhi.init(0,maxPt,mPt, -maxPsi,maxPsi,nPsi, "");
  his_onePtPhi.init(0,maxPt,mPt, -maxPsi,maxPsi,nPsi,
                  "../out/track/sim_effOne_ptPhi_"+partName+".his"); 

  //
  theResolutions = new Resolutions(1,partName); // for simu, pt resolution

  //
  string path = "../simu/simulation_"+partName+".dat.gz";
  cerr << Helper::col(2) << " reading from " << path << Helper::col();
  igzstream file(path.c_str());

  char code;
  vector<Part> genParts, simParts, recParts;

  int clusterCountBpix;
  unsigned int layerCodes;
  int nBpixTracks, nPixelTracks;

  while(file >> code)
  {
    if(code == 'g' && genParts.size() > 0) // evaluate event
    {
      if(genParts.size() != 1) exit(1);

      const Part & gen = genParts[0];

      his_genTrack.fill({gen.eta, gen.pt});  

      // is well reconstructed?
      bool isRecod = (recParts.size() == 1);

      if(isRecod)
        his_oneTrack.fill({gen.eta, gen.pt});  
      
      if(recParts.size() > 1)
        his_mulTrack.fill({gen.eta, gen.pt});  

      //
      int layerCountBpix = getLayerCountBpix(layerCodes);

      bool wouldFireHlt =
           ( (layerCountBpix >= minLayers && clusterCountBpix >= minClusters) ||
              nBpixTracks  >= minBpixTracks ||
              nPixelTracks >= minPixelTracks);   

      // vs phi
      if(wouldFireHlt)
      {
        if(gen.pt < 1.) // pT < 1 GeV
        {
          his_genEtaPhi.fill({gen.eta, gen.phi});

          if(isRecod)
          his_oneEtaPhi.fill({gen.eta, gen.phi});
        }

        his_genPtPhi.fill({gen.pt, gen.phi});

        if(isRecod)
        his_onePtPhi.fill({gen.pt, gen.phi});
      }

      // exactly one reconstructed
      if(isRecod)
      if(wouldFireHlt)
        his_hltTrack.fill({gen.eta, gen.pt});

      // store info
      vector<short int> index;
      if(getIndex(gen, index))
      {
        SimTrackInfo & info = simInfo[index];
    
                           info.nall++; 
        if(isRecod)
        {
                           info.nrec++; 

          if(wouldFireHlt) info.nhlt++; 
          else
          {
            // here nBpixTracks=0, nPixelTracks=0

            // set to min
            if(clusterCountBpix >= minClusters)
               clusterCountBpix  = minClusters;

            pair<unsigned char, unsigned char> code(
                 layerCodes, clusterCountBpix);

            info.noHltCodes[code]++;
          }
        }
      }

      // for pt resolution
      if(isRecod)
      {
        const Part & rec = recParts[0];

        theResolutions->collectForPt({gen.eta, gen.pt,
                                      rec.pt - gen.pt});
      }

      // clear
      genParts.clear(); simParts.clear(); recParts.clear();

      Helper::propeller(nEvents++,1e+6);
    }

    if(code == 'g') // genparticle
    {
      Part part;
      file >> part.pid >> part.eta >> part.pt >> part.phi;
      genParts.push_back(part);
    }

    if(code == 'h') // information for the HLT trigger
    {
      file >> clusterCountBpix >> layerCodes >> nBpixTracks >> nPixelTracks;
    }

    if(code == 't') // high-level pixel-cluster and pixel-track trigger
    {
      int thisWouldFireHlt;
      file >> thisWouldFireHlt;
    }

    if(code == 'r') // rectrack
    {
      Part part;
      file >> part.pid >> part.eta >> part.pt >> part.phi; // here pid = q

      if(isCompatible(genParts[0], part))
        recParts.push_back(part);
    }

    if(code == 's') // simtrack
    {
      Part part;
      file >> part.pid >> part.eta >> part.pt >> part.phi
              >> part.vz >> part.nMatched;

      if(abs(part.pid) != 11 &&
         abs(part.pid) != 13) // not secondary e, mu
        simParts.push_back(part);
    }
  }

  file.close();

  cerr << " [read " << nEvents << " events]" << endl;

  //
  his_hltTrack.div(his_oneTrack);

  his_oneTrack.div(his_genTrack);
  his_mulTrack.div(his_genTrack);

  his_oneEtaPhi.div(his_genEtaPhi);
  his_onePtPhi.div(his_genPtPhi);

  // write
  char name[256];
  sprintf(name,"../out/track/simTable_%s.out.gz",partName.c_str());
  ogzstream fileOut(name);

  sprintf(name,"../out/track/simTable_%s.zeros",partName.c_str());
  ofstream fileZeros(name);

  for(short int ieta = 0; ieta < nEta; ieta++)
  for(short int ipt  = 0; ipt  < mPt ; ipt++)
  for(short int ipsi = 0; ipsi < nPsi; ipsi++)
  {
    vector<short int> index = {ieta,ipt,ipsi};
    SimTrackInfo & info = simInfo[index];

    fileOut << " " << info.nall << " " << info.nrec << " " << info.nhlt;

    if(info.nrec < info.nall*0.1)
      fileZeros << " " << -maxEta + (ieta + 0.5)/nEta * 2*maxEta
                << " " <<           (ipt  + 0.5)/mPt  *   maxPt
                << " " << -maxPsi + (ipsi + 0.5)/nPsi * 2*maxPsi
                << endl;

    for(auto & c : info.noHltCodes)
      fileOut << "  " << int(c.first.first)  // code (layerCodes)
               << " " << int(c.first.second) // code (clusterCountBpix)
               << " " <<     c.second;       // times

    fileOut << endl;
  }

  fileOut.close();

  fileZeros.close();

  //
  delete theResolutions;
}

/*****************************************************************************/
void TrkEffic::readSimTable(const string & partName,
                            map<vector<short int>,SimTrackInfo> & simInfo)
{
  // read
  char name[256];
  sprintf(name,"../out/track/simTable_%s.out.gz",partName.c_str());
  igzstream file(name);

  cerr << ".";

  int nRead = 0;

  //
  for(short int ieta = 0; ieta < nEta; ieta++)
  for(short int ipt  = 0; ipt  < mPt ; ipt++)
  for(short int ipsi = 0; ipsi < nPsi; ipsi++)
  {
    vector<short int> index = {ieta,ipt,ipsi};
    SimTrackInfo & info = simInfo[index];

    file >> info.nall >> info.nrec >> info.nhlt;

    if(!file.eof()) nRead++;

    int num = info.nrec - info.nhlt;

    int sum = 0;

    while(sum < num)
    {
      int layerCodes, clusterCountBpix, times;
      file >> layerCodes >> clusterCountBpix >> times;

      pair<unsigned int,unsigned int> code(layerCodes, clusterCountBpix);

      info.noHltCodes[code] = times;

      sum += times;
    }
  }

  file.close();
}

/*****************************************************************************/
void TrkEffic::readSimTables()
{
  cerr << Helper::col(2) << " reading HLT+tracking effic" << Helper::col();

  readSimTable("pip", simInfo[pion][0]);
  readSimTable("pim", simInfo[pion][1]);

  readSimTable("kap", simInfo[kaon][0]);
  readSimTable("kam", simInfo[kaon][1]);

  readSimTable("prp", simInfo[prot][0]);
  readSimTable("prm", simInfo[prot][1]);

  cerr << " [done]" << endl;
}

/*****************************************************************************/
double TrkEffic::getEfficiency(Part p[2], int type)
{
  vector<short int> index[2];

  if(getIndex(p[0],index[0]) &&
     getIndex(p[1],index[1]))
  {
    SimTrackInfo & pos = simInfo[type][0][index[0]];
    SimTrackInfo & neg = simInfo[type][1][index[1]];

    int sum = 0;
    int num = pos.nall * neg.nall;

    // rec-nohlt and rec+hlt
    sum += pos.nhlt * neg.nhlt;

    // rec-nohlt and rec+hlt
    sum += (pos.nrec - pos.nhlt) * (neg.nhlt);

    // rec+hlt and rec-nohlt
    sum += (pos.nhlt) * (neg.nrec - neg.nhlt);

    // rec-nohlt and rec-nohlt [both reconstructed]
    for(auto & mpos : pos.noHltCodes)
    for(auto & mneg : neg.noHltCodes)
    {
      const pair<unsigned char,unsigned char> & cpos = mpos.first;
      const pair<unsigned char,unsigned char> & cneg = mneg.first;

      const int & npos = mpos.second;
      const int & nneg = mneg.second;

      int layerCodes       = cpos.first  | cneg.first;
      int clusterCountBpix = cpos.second + cpos.second;

      int layerCountBpix = getLayerCountBpix(layerCodes);

      bool wouldFireHlt =
        (layerCountBpix >= minLayers && clusterCountBpix >= minClusters);

      if(wouldFireHlt) sum += npos * nneg;
   }

    return double(sum) / num;
  }
  else
    return 0;
}

/*****************************************************************************/
double TrkEffic::getEfficiency(const vector<CmTrack> & tracks, int type)
{
  if(tracks.size() == 2)
  {
    Part parts[2];

    for(int i = 0; i < 2; i++)
    {
      parts[i].eta =  tracks[i].p.eta();
      parts[i].pt  = (tracks[i].p.trans() < 2 ?
                      tracks[i].p.trans() : 2 - 1e-3); // - small
      parts[i].phi =  tracks[i].p.phi();
    }

    return getEfficiency(parts, type);
  }
  else
    exit(1);
}

/*****************************************************************************/
double TrkEffic::getEfficiency(const Vector3 & p3,
                               const Vector3 & p4, int type)
{
  Part parts[2];

  parts[0].eta =  p3.eta();
  parts[0].pt  = (p3.trans() < 2 ?
                  p3.trans() : 2 - 1e-3); // - small
  parts[0].phi =  p3.phi();

  parts[1].eta =  p4.eta();
  parts[1].pt  = (p4.trans() < 2 ?
                  p4.trans() : 2 - 1e-3); // - small
  parts[1].phi =  p4.phi();

  return getEfficiency(parts, type);
}

