#include "../interface/ProcessData.h"

#include <thread>
#include <algorithm>
#include <experimental/filesystem>

#include "../interface/RpPrint.h"
#include "../interface/RpPat.h"
#include "../interface/RpFit.h"
#include "../interface/RpEffic.h"
#include "../interface/RpReco.h"
#include "../interface/RpVeto.h"

#include "../interface/ParticleId.h"
#include "../interface/TrkEffic.h"
#include "../interface/TrkReco.h"
#include "../interface/Selection.h"

#include "../interface/Physics.h"

#include "../interface/Helper.h"
#include "../interface/gzstream.h"

#define sqr(x) ((x)*(x))

using namespace std;

#undef Debug

/*****************************************************************************/
ProcessData::ProcessData(
  const string & data, int nCmTra, int nEvery, int nCpus,
                       const vector<string> & flags) :
  data(data), nCmTra(nCmTra), nEvery(nEvery), nCpus(nCpus), flags(flags)
{
  //
  cerr << Helper::col(6) << " [using " << data << "]" << Helper::col() << endl;

  //
  if(has("-rpPats" ))
  {
    theRpPat   = new RpPat(0, nCmTra);
    cerr << Helper::col(1)
         << " collecting roman pot strip patterns.."
         << Helper::col() << endl;
    return;
  }

  //
  if(has("-rpEffic"))
  {
    theRpEffic = new RpEffic(0, nCmTra);
    cerr << Helper::col(1)
         << " collecting roman pot strip efficiencies.."
         << Helper::col() << endl;
    return;
  }

  if(has("-rpVeto"))
  {
    theRpReco     = new RpReco(nCmTra);
    theRpVeto     = new RpVeto(1, nCmTra); // collect
    return;
  }

  //
  if(has("-rpReco") || has("-doPhys"))
  {
    theRpReco     = new RpReco(nCmTra);
    theRpVeto     = new RpVeto(-1, nCmTra); // do nothing (elasticMask)

    if(nCmTra == 2)
    {
      theParticleId = new ParticleId(1);      // collect relsig
      theTrkEffic   = new TrkEffic(1);        // read simTables for getEffic
      theTrkReco    = new TrkReco(1);         // collect histos
    }

    theSelection  = new Selection(nCmTra);

    cerr << Helper::col(1)
         << " reconstructing roman pot tracklets.."
         << Helper::col() << endl;
  }

  if(has("-doPhys"))
  {
    thePhysics = new Physics(nCmTra);

    cerr << Helper::col(1)
         << " doing physics.."
         << Helper::col() << endl;
  }
}

/*****************************************************************************/
ProcessData::~ProcessData()
{
  if(has("-rpPats" )) delete theRpPat;   // write strip patterns
  if(has("-rpEffic")) delete theRpEffic; // write strip effic
  if(has("-rpVeto" )) delete theRpVeto;  // write !veto (py1,py2) efficiency
  if(has("-rpReco" ) || has("-doPhys" )) // FIXME
  {
    delete theRpReco;  // write histos

    if(nCmTra == 2)
    {
      delete theParticleId;
      delete theTrkReco;
      delete theSelection;
    }

    //
    if(nCmTra == 2)
    {
      cerr << Helper::col(5) << " writing taken events.." << Helper::col();

      ofstream file("../out/lumi/takenEvents.dat");
      int orun = -1;

      for(auto & a : takenEvents)
      {
        const int & run = a.first.first;

        if(run != orun)
        {
          if(orun != -1) file << endl << endl;
          orun = run;
        }

        file << " " << a.first.first
             << " " << a.first.second
             << " ";
  
        for(int t = 0; t < nTopos; t++)
          file << " " << a.second[t];

        file << endl;
      }
      file.close();

      cerr << " [done]" << endl;
    }
  }

  if(has("-doPhys" )) delete thePhysics;
}

/*****************************************************************************/
bool ProcessData::has(const string & flag)
{
  return (find(flags.begin(), flags.end(), flag) != flags.end());
}

/*****************************************************************************/
bool ProcessData::readEvent(igzstream & fileIn, int nCmTra, Event & event)
{
  // clear
  event.run = -1;
  event.ls  = -1;
  event.bx  = -1;

  event.rpDets.clear();
  event.rpTracks.clear();
  event.prTracks.clear();
  event.cmTracks.clear();

  const bool rpPrint = has("-rpPrint");

  // read
  if(fileIn >> event.run >> event.ls >> event.bx)
  {
    int nDets;
    if(rpPrint)
      fileIn >> nDets;

    int nRpTracks;
    fileIn >> nRpTracks;

    if(rpPrint)
    {
      // roman pot dets
      for(int i = 0; i < nDets; i++)
      {
        RpDet det;
        int nClus;
        fileIn >> det.arm >> det.sta >> det.rpt >> det.pla >> nClus;

        det.sta /= 2; // 0 2 -> 0 1
        det.rpt %= 4; // 4 5 -> 0 1

        det.uv = det.pla % 2;
        det.pla /= 2;

        for(int j = 0; j < nClus; j++)
        {
          double strip;
          int width;

          fileIn >> strip >> width;
          if(strip == -1) strip = empty; // -1 -> -99
          det.clus.push_back(strip);
        }

        event.rpDets.push_back(det);
      }
    }

    // roman pot tracks
    for(int i = 0; i < nRpTracks; i++)
    {
      RpTrack track;
      RpDet   & det = track.det;
      Vector2 & pos = track.pos;

      fileIn >> det.arm >> det.sta >> det.rpt;

      bool ok = (det.sta == 0 || det.sta == 2) ||
                (det.rpt == 4 || det.rpt == 5);

      det.sta /= 2; // 0 2 -> 0 1
      det.rpt %= 4; // 4 5 -> 0 1

      float sigx, sigy;
      fileIn >> pos.x; if(!rpPrint) fileIn >> sigx;
      fileIn >> pos.y; if(!rpPrint) fileIn >> sigy;

      pos.cxx = sqr(sigx);
      pos.cyy = sqr(sigy);

      for(short int pla = 0; pla < nPlanes; pla++)
      for(short int uv = 0; uv < 2; uv++)
      {
        double strip;
        fileIn >> strip;

        if(det.print() == "0|0|0" && uv == 0 && pla == 2) // FIXME VERY
          strip = -1;

        if(strip == -1) strip = empty; // -1 -> -99
        track.clus[uv][pla] = strip;
      }

      //
      if(ok) event.rpTracks.push_back(track);
    }

    // cms track vertex
    if(!rpPrint)
    if(nCmTra > 0 && nCmTra != 9)
    {
      CmVertex & vtx = event.cmVertex;
      fileIn >> vtx.x.val
             >> vtx.y.val
             >> vtx.z.val;

      float cxx,cyy,czz;
      fileIn >> cxx >> cyy >> czz;
      vtx.x.sig = sqrt(cxx);
      vtx.y.sig = sqrt(cyy);
      vtx.z.sig = sqrt(czz);

      fileIn >> vtx.chi2 >> vtx.ndf;
    }

    // cms tracks
    if(nCmTra > 0 && nCmTra != 9)
    for(int i = 0; i < nCmTra; i++)
    {
      CmTrack track;
      string s;
      int nhits;

      fileIn
         >> track.eps >> s >> nhits
         >> track.q
         >> track.p.x >> track.p.y >> track.p.z;

      fileIn >> track.dt.val; if(!rpPrint) fileIn >> track.dt.sig;
      fileIn >> track.dz.val; if(!rpPrint) fileIn >> track.dz.sig;
      fileIn >> track.chi2 >> track.ndf;

      if(!rpPrint) fileIn >> track.sigpt;

      if(!fileIn.eof())
      {
        if(s != "inf" && s != "nan") track.sigma = stof(s);
                                else track.sigma = 1/0.;
        event.cmTracks.push_back(track);
      }
      else 
        return false;
    }

    // swap to +, -
    if(event.cmTracks.size() == 2)
    if(event.cmTracks[0].q < 0 && event.cmTracks[1].q > 0)
      swap(event.cmTracks[0], event.cmTracks[1]);

    return true;
  }
  else
    return false;
}

/*****************************************************************************/
bool ProcessData::processEvent(Event & event)
{
  if(has("-rpPrint"))
  { 
    RpPrint theRpPrint;
    theRpPrint.print(event.rpDets, event.rpTracks);

    getchar();

    return true;
  }

#ifdef Debug
cerr << " --------------------" << endl;
cerr << " a0 process" << endl;
#endif

  // need exactly four tracklets, proper topology
  if(theRpReco->getTopo(event.rpTracks, event.topo))
  {
#ifdef Debug
cerr << " a1 topo " << event.topo << endl;
#endif

    // collect patterns
    if(has("-rpPats"))
    {
      theRpPat->collectPatterns(event.rpTracks);
      return true;
    }

    // collect efficiency
    if(has("-rpEffic"))
    {
      theRpEffic->collectStrip(event.rpTracks, event.run);
      return true;
    }

    if(has("-rpVeto"))
    if(theRpReco->reconstruct(event.topo, event))
    {
      theRpVeto->combineTracks(event.rpTracks, event.prTracks, event.topo);
      return true;
    }

    if(has("-rpReco") || has("-doPhys"))
    {
      // roman pots reconstruction | problems with global vars in RpFit (my_*)
      // get event.rpWeight!
      // if hit locations |y| are not within limits, skip event
      if(!theRpReco->reconstruct(event.topo, event))
        return false;

#ifdef Debug
cerr << " a2" << endl;
#endif

      //////////////////////////////////////////////////////////////
      // event filters
      bool isAlright = true;
      int problem = -1;
      event.type = unknown;

      mtx_his.lock(); // lock
      if(nCmTra == 2)
      {
        theTrkReco->collectPars(event.cmTracks); // collect

        // identify pair
        event.type = theParticleId->identifyPair(event.cmTracks);
        if(event.type == unknown)
        { isAlright = false; problem = 1; }

        // guess pid, looper removal
        int typeGuess = (event.type != unknown ? event.type : pion);

        if(theSelection->areCmTracksLooper(event,typeGuess, event.rpWeight))
        { isAlright = false; problem = 2; }

        // primary central hadros
        if(!theSelection->areCmTracksPrimaries(event, event.rpWeight))
        { isAlright = false; problem = 3; }
      }

      // primary scattered protons
      if(!theSelection->areRpProtonsPrimaries(event, event.rpWeight))
        { isAlright = false; problem = 4; }

#ifdef Debug
cerr << " a3 alright = " << isAlright << " " << problem << endl;
#endif

      //////////////////////////////////////////////////////////////
      // filtered events only
      // central hadrons are (not looper, compat with IP, identified, same type)
      // scattrd protons are (compat with each other at IP)
      if(isAlright)
      {
        // classify event
        event.cat = theSelection->classifyEvent(event, event.topo,
                                                event.rpWeight);

        event.sign = 0;
        if(event.cat == signal)   event.sign =  1;
        if(event.cat == sideband) event.sign = -1;

        if(nCmTra == 2)
        {
          // collect
          theParticleId->processEnergyLoss(event.cmTracks,
                                           event.cat, event.type);

#ifdef Debug
cerr << " a4 sign = " << event.sign << endl;
#endif
         
          // signal or sideband
          if(event.sign != 0)
          {
            theTrkReco->collectChi2(event.cmTracks, event.type); // collect

#ifdef Debug
cerr << " a5 good effic   = " << theTrkReco->hasGoodEffic(event.cmTracks) << endl;
cerr << " a6 elastic mask = " << theRpVeto->elasticMask(event.prTracks, event.topo)<< endl;
#endif

            // larger than real mask ok (p1y,p2y) for takenEvents in TB,BT
            if(theTrkReco->hasGoodEffic(event.cmTracks))
            if(!theRpVeto->elasticMask(event.prTracks, event.topo))
            { // collect for luminosity comparison
              pair<int,int> loc(event.run,event.ls);

#ifdef Debug
cerr << " a7 rpWeight = " << event.rpWeight << endl;
#endif

              // get pair efficiency
              double trkEff =
                theTrkEffic->getEfficiency(event.cmTracks, event.type);

              if(trkEff < 0.05) trkEff = 0.05; // FIXME PARAMETER

              // calculate track weight
              double trkWeight = 1/trkEff;

              // sum
              takenEvents[loc][event.topo] +=
                 event.sign * event.rpWeight * trkWeight;
            }

            // do physics
            if(has("-doPhys"))
              thePhysics->process(event);

#ifdef Debug
getchar();
#endif
          }
        }
      }

      mtx_his.unlock(); // unlock

      return true;
    }
  }

  return false; 
}

/*****************************************************************************/
void ProcessData::readStreams(const vector<string> & dataStreams)
{
  // take all streams
  for(auto & dataStream : dataStreams)
  {
    const string path = "../data/"+data+"/"+dataStream+"/";

    if(experimental::filesystem::exists(path))
    {
      int nFiles = 0;

      mtx_cerr.lock();
      cerr << " " << dataStream;
      mtx_cerr.unlock();

      // all files
      for(const auto & entry :
          experimental::filesystem::directory_iterator(path))
      if(nFiles++ % nEvery == 0)
      {
        igzstream fileIn(entry.path().string().c_str());

        //
        Event event;

        while(readEvent(fileIn, nCmTra, event))
        {
//          if(event.run == 319311)
          if(processEvent(event))
          {
            mtx_count.lock();
            nEvents++; Helper::propeller(nEvents, 100000);
            mtx_count.unlock();
          }
        }

        fileIn.close();
      }
    }
  }
}

/*****************************************************************************/
void ProcessData::readData()
{
  nEvents = 0;

  // setup lists
  vector<vector<string>> lists =
    Helper::getStreams(nCpus, has("-rpVeto")); // paraOnly?

  cerr << " reading";

  // initialize threads
  vector<thread> threads;
  for(auto & list : lists)
    threads.push_back(thread(&ProcessData::readStreams,this,list));

  // start and wait
  for(auto & thread : threads)
    thread.join();

  cerr << endl;

  // final report
  if(nEvents > 1e+6)
    cerr << " [read " << int(nEvents/1e+6 * 100)/100. << " M events]";
  else
    cerr << " [read " <<     nEvents                  <<   " events]";

  cerr << endl;
}

