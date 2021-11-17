#include "../interface/RpSimu.h"

#include "../interface/Parameters.h"

#include "../interface/RpFit.h"
#include "../interface/RpEffic.h"
#include "../interface/Random.h"
#include "../interface/Histo.h"
#include "../interface/Helper.h"

#include "../interface/gzstream.h"

#include <thread>

using namespace std;

/*****************************************************************************/
RpSimu::RpSimu(int nCpus, int nEvents) : nCpus(nCpus), nEvents(nEvents)
{
  int nCmTra = 2;

  const map<RpDet, map<vector<double>,int>> dummy;
  theRpFit   = new RpFit(dummy,nCmTra,1); // read shifts
  theRpEffic = new RpEffic(1,  nCmTra);   // read strip efficiencies
  theRandom  = new Random();
}

/*****************************************************************************/
RpSimu::~RpSimu()
{
  delete theRpFit;
  delete theRpEffic;
  delete theRandom;
}

/*****************************************************************************/
bool RpSimu::isSingle(double x)
{
  return (x == int(x));
}

/*****************************************************************************/
double RpSimu::generateHit()
{
/*
  // double gaussian
  const double m = 3*(nStrips/4);

  double x,p;

  do
  {
    x = theRandom->getFlat() * nStrips;

    double s = (x < m ? 100 : 50);

    double q = (x-m)/s;

    p = exp(-q*q/2);
  }
  while(theRandom->getFlat() > p);
*/

  // flat
  mtx_rnd.lock();
  double x = theRandom->getFlat() * nStrips;
  mtx_rnd.unlock();

  return x;
}

/*****************************************************************************/
bool RpSimu::generateRpTrack(int run, RpTrack & track)
{
  vector<double> a,b;

  for(int uv = 0; uv < 2; uv++)
  {
    mtx_rnd.lock();
    a.push_back(theRandom->getFlat(-maxSlope,maxSlope)); // slope PARAMETER
    mtx_rnd.unlock();

    b.push_back(generateHit());                // intercept
  }

  RpDet det = track.det;

  // choose a layer for nuclear interaction
  int sLay = 99;
  double sNuc = 0;

  mtx_rnd.lock();
  if(theRandom->getFlat() < 0.20) // PARAMETER
  {
    sLay = int(theRandom->getFlat(0,10)); // position of scatter
    sNuc = 0.7*theRandom->getGauss(); // shift per layer distance PARAMETER
  }
  mtx_rnd.unlock();

  // generate tracklets
  for(int uv  = 0; uv  < 2;       uv++ )
  for(int pla = 0; pla < nPlanes; pla++)
  {
    det.uv = uv;

    // shift due to nuclear collision in layer lNucl, simple model both u/v
    int dLay = 2*pla+uv - sLay;
    double nuclShift = (dLay > 0 ? dLay * sNuc : 0);

    // get measured (single or double)
    double strip = theRpFit->getStrip(det,pla,a[uv],b[uv] + nuclShift);

    // use strip efficiency
    if(strip >= 0 && strip <= nStrips-1)
    {
      if(isSingle(strip))
      { // single
        int istrip = round(strip);

        det.pla = pla;
        det.str = istrip; double effic = theRpEffic->getStripEffic(run, det);
        det.pla = -1; det.str = -1;

        mtx_rnd.lock();
        bool isEff = (theRandom->getFlat() < effic);
        mtx_rnd.unlock();

        if(!isEff) strip = -1; // if inefficienct;
      }
      else
      { // double
        int istrip = int(strip);

        det.pla = pla;
        vector<double> effic(2);
        det.str = istrip;   effic[0] = theRpEffic->getStripEffic(run, det);
        det.str = istrip+1; effic[1] = theRpEffic->getStripEffic(run, det);
        det.pla = -1; det.str = -1;

        vector<bool> isEff(2);
        mtx_rnd.lock();
        for(int i = 0; i < 2; i++)
          isEff[i] = (theRandom->getFlat() < effic[i]); 
        mtx_rnd.unlock();

        if( isEff[0] && !isEff[1]) strip = istrip;
        if(!isEff[0] &&  isEff[1]) strip = istrip+1;
        if(!isEff[0] && !isEff[1]) strip = -1;
      }
    }
    else
    {
      strip = -1;
    }

    track.clus[uv][pla] = strip;
  }

  // check if we have enough recorded hits per orientation
  bool ok = true;

  for(int uv = 0; uv < 2; uv++)
  {
    det.uv = uv;
    int nhits = 0;

    // count hits
    map<int,int> ntrig;
    const int tStrips = 16; // numer of "trigger strips"
    const int wTrig = nStrips / tStrips;

    for(int pla = 0; pla < nPlanes; pla++)
    {
      const double & str = track.clus[uv][pla];

      if(str != -1)
      {
        nhits++;

        double red = str - int(str / wTrig) * wTrig;

        if(red < wTrig-0.5) // must not be 31.5
          ntrig[int(str) / wTrig]++;
      }
    }

    bool trig_ok = false;
    for(int i = 0; i < tStrips; i++)
    if(ntrig[i] >= 3) trig_ok = true;

    // enough hits and for far pots enough hits in the same "trigger strip"
    bool isDetected = (nhits >= 3 && (det.sta == 0 || trig_ok)); // 

    if(!isDetected) ok = false; // need at least three hits

    if(!toFile)
    {
      mtx_his.lock();
      if(isDetected)
        his_effic[det].fill({b[uv],a[uv]});

      den_effic[det].fill({b[uv],a[uv]});
      mtx_his.unlock();
    }
  }

  return ok;
}

/*****************************************************************************/
void RpSimu::printRpTrack(ogzstream & file, RpTrack & track)
{
  const RpDet & det = track.det;

  file << " " << det.arm << " " << (det.sta << 1) << " " << 4 + det.rpt
       << " 0. 0. 0. 0.";

  for(int pla = 0; pla < nPlanes; pla++)
  for(int uv = 0; uv < 2; uv++)
    file << " " << track.clus[uv][pla];

  file << endl;
}

/*****************************************************************************/
void RpSimu::generateStream(int run, const vector<string> & dataStreams)
{
  int mEvents = (nEvents * dataStreams.size()) / 8;

  for(auto & dataStream : dataStreams)
  {
    ogzstream file;

    if(toFile)
    {
      char name[256];
      sprintf(name, "../data/9part/%s/%d.dat.gz",dataStream.c_str(),run);
      file.open(name);
    }

    for(int event = 0; event < mEvents; event++)
    {
      mtx_rnd.lock();
      int topo = (theRandom->getFlat() < 0.5 ? 0 : 1); // TB, BT
      mtx_rnd.unlock();

      //
      int nok_far = 0;

      vector<RpTrack> rpTracks;

      for(int arm = 0; arm < 2; arm++)
      {
        int rpt = (arm ^ topo); // TB a1-> T, a2 -> B | BT a1 -> B, a2 -> T

        for(int sta = 0; sta < 2; sta++)
        {
          RpDet det;
          det.arm = arm; det.rpt = rpt; det.sta = sta;

          RpTrack track;
          track.det = det;

          track.pos.x = 0;
          track.pos.y = 0;

          track.pos.cxx = 0;
          track.pos.cyy = 0;
          track.pos.cxy = 0;

          if(generateRpTrack(run, track))
          {
            rpTracks.push_back(track);

            if(sta == 1) nok_far++; // count far
          }
        }
      }

      if(nok_far == 2) // L1 triggered
      if(toFile)
      {
        // event header
        const int ls = 1;
        const int bx = 2;

        file << run << " " << ls << " " << bx << " " << rpTracks.size() << endl;

        for(auto & track : rpTracks)
          printRpTrack(file,track);
      }

      mtx_count.lock();
      iEvents++; Helper::propeller(iEvents, nEvents/10);
      mtx_count.unlock();
    }

    if(toFile)
      file.close();
  }
}

/*****************************************************************************/
void RpSimu::generateData(int nCmTra, int run)
{
  toFile = (nCmTra == 9);

//  cerr << (toFile ? " generating to file.."
//                  : " filling efficiency histos..") << endl;

  if(!toFile)
  {
    // init histos
    his_effic.clear();
    den_effic.clear();

    for(int a = 0; a < 2; a++)
    for(int s = 0; s < 2; s++)
    for(int r = 0; r < 2; r++)
    for(int u = 0; u < 2; u++)
    {
      RpDet det;
      det.arm = a;
      det.sta = s;
      det.rpt = r;
      det.uv  = u;

      char name[256];
      sprintf(name,"../out/rp/%dpart/groupEffic/%d_%d%d%d_%d.his",
                   nCmTra, run,a,s,r, u);

      his_effic[det].init(0,nStrips,nStrips, -maxSlope,maxSlope,binSlope, name);
      den_effic[det].init(0,nStrips,nStrips, -maxSlope,maxSlope,binSlope, "");
    }
  }

  //
  cerr << "  run " << run;

  iEvents = 0;

  // lists for threads
  vector<vector<string>> lists =
    Helper::getStreams(nCpus, false);

  // initialize and start threads
  vector<thread> threads;
  for(auto & list : lists)
    threads.push_back(thread(&RpSimu::generateStream,this, run,list));

  for(auto & thread : threads)
    thread.join();

  //
  if(!toFile)
  {
    cerr << Helper::col(6) << " writing.."
         << Helper::col();

    for(auto & m : his_effic)
    {
      const RpDet & det = m.first;

      his_effic[det].div(den_effic[det]);
      his_effic[det].write();
    }
  }

  //
  cerr << " [done]" << endl;
}

