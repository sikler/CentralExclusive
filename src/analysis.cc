#include "../interface/Helper.h"
#include "../interface/RpPat.h"
#include "../interface/RpFit.h"
#include "../interface/RpSimu.h"
#include "../interface/RpVeto.h"
#include "../interface/ProcessData.h"
#include "../interface/Luminosity.h"
#include "../interface/TrkEffic.h"
#include "../interface/TrkReco.h"
#include "../interface/ParticleId.h"

#include <iostream>
#include <cstring>

using namespace std;

int nCmTra,nEvery,nCpus;
string data;
vector<string> flags;

/*****************************************************************************/
void options(int arg, char *arc[])
{
  nCmTra = -1; // initialization needed
  nEvery =  1; // default
  nCpus  =  1; // default

  for(int i = 1; i < arg; i++)
  {
    if(strcmp(arc[i],"-nCmTra") == 0) { nCmTra = atoi(arc[++i]); continue; }
    if(strcmp(arc[i],"-nEvery") == 0) { nEvery = atoi(arc[++i]); continue; }
    if(strcmp(arc[i],"-nCpus")  == 0) { nCpus  = atoi(arc[++i]); continue; }

    flags.push_back(arc[i]); // collect other flags
  }

  //
  if(nCmTra != -1)
    data = to_string(nCmTra) + "part";
  else 
  {
    cerr << " error: need nCmTra != -1" << endl;
    exit(1);
  }
}

/*****************************************************************************/
int main(int arg, char *arc[])
{
  options(arg, arc);

  /////////////////////////////////////////////////////////////////////////////
  // -data
  if(strcmp(arc[1],"-data") == 0)
  {
    if(strcmp(arc[2],"-rpPrint") == 0)
      data = "2print";

    ProcessData theProcessData(data, nCmTra,nEvery,nCpus, flags);
    theProcessData.readData();
  }

  /////////////////////////////////////////////////////////////////////////////
  // -lumi
  if(strcmp(arc[1],"-lumi") == 0)
  {
    if(strcmp(arc[2],"-calcWeights") == 0)
    {
      cerr << Helper::col(1)
           << " luminosity: calculate weights"
           << Helper::col() << endl;

      Luminosity theLuminosity(true ,false); // readBunchInfo, readDataInfo
      theLuminosity.calcWeights();
    }

    if(strcmp(arc[2],"-postProcess") == 0)
    {
      cerr << Helper::col(1)
           << " luminosity: post-process"
           << Helper::col() << endl;

      Luminosity theLuminosity(true ,true ); // readBunchInfo, readDataInfo
      theLuminosity.postProcess();
    }
  }

  /////////////////////////////////////////////////////////////////////////////
  // -calc
  if(strcmp(arc[1],"-calc") == 0)
  {
    // fit roman pot shifts
    if(strcmp(arc[2],"-rpShifts") == 0)
    {
      RpPat theRpPat(2, nCmTra); // read dominant (>= 100)
      const auto & patterns = theRpPat.getPatterns(false);

      //
      RpFit theRpFit(patterns,nCmTra, 0); // nothing

      theRpFit.optimizeShifts(true);  // collect chi2
      theRpFit.optimizeShifts(false); // count zeros
    }

    // fit roman pot patterns
    if(strcmp(arc[2],"-rpPats") == 0)
    {
      RpPat theRpPat(3, nCmTra); // no 1nTu3 hit
      const auto & patterns = theRpPat.getPatterns(true); // add holes

      //
      RpFit theRpFit(patterns,nCmTra, 1); // read shifts
      theRpFit.fitPatterns();
    }

    // generate roman pot data (nCmTra = 9) to file
    if(strcmp(arc[2],"-rpGene") == 0)
    {
      const int nEvents = 1e+6; // PARAMETER

      cerr << Helper::col(1)
           << " generate roman pot data to file.."
           << Helper::col()
           << " ("<<nEvents/nEvery/1e+6 << "M events per data stream)"
           << endl;

      RpSimu theRpSimu(nCpus,nEvents/nEvery);

      const int run = 319311; // efficiencies taken from this run
      theRpSimu.generateData(nCmTra, run);
    }

    // calculate roman pot tracklet (strip-group) efficiencies
    if(strcmp(arc[2],"-rpGroup") == 0)
    {
      const int nEvents = 2e+6; // PARAMETER

      cerr << Helper::col(1)
           << " calculate roman pot tracklet (strip-group) efficiencies.."
           << Helper::col()
           << " ("<<nEvents/nEvery/1e+6 << "M events per run)"
           << endl;

      RpSimu theRpSimu(nCpus,nEvents/nEvery);

      for(auto & run : runs)
        theRpSimu.generateData(nCmTra, run);
    }

    // calculate roman pot angular efficiency
    if(strcmp(arc[2],"-rpAngEff") == 0)
    {
      RpVeto theRpVeto(0, nCmTra); // read, calculate angular efficiency
      theRpVeto.calculateAngularEfficiency();
    }

    // calculate combined HLT/tracking/PID efficiency
    if(strcmp(arc[2],"-cmEff") == 0)
    {
      int type = unknown;
      for(int i = 0; i < nParts; i++)
        if(arc[3] == parts[i]) type = i;

      TrkReco theTrkReco(0);
      theTrkReco.prepareCmTables(type);
    }
  }

  /////////////////////////////////////////////////////////////////////////////
  // -pidDemo
  if(strcmp(arc[1],"-pidDemo") == 0)
  {
    cerr << Helper::col(1) << " preparing PID demo for " << arc[2]
         << Helper::col() <<endl;

    ParticleId theParticleId(2); // use relSigma

    int type = unknown;
    for(int i = 0; i < nParts; i++)
      if(arc[2] == parts[i]) type = i;

    theParticleId.prepareDemo(type);
  }

  /////////////////////////////////////////////////////////////////////////////
  // -trkEff
  // Tracking efficiency
  if(strcmp(arc[1],"-trkEff") == 0)
  {
    TrkEffic theTrkEffic(0); // process simulation
    theTrkEffic.processSimulation(arc[2]);
  }

}

