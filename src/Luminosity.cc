#include "../interface/Luminosity.h"

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>

#include "../interface/Parameters.h"
#include "../interface/Helper.h"
#include "../interface/gzstream.h"

using namespace std;

enum { filled = 1, selected = 2 };

const double sig_vis  = 79e-3;   // [b] PARAMETER

const double sig_diag = 15.2e-6; // [b] PARAMETER
const double sig_para = 12.4e-6; // [b] PARAMETER

/*****************************************************************************/
Luminosity::Luminosity() // read effective luminosity weights
{
}

/*****************************************************************************/
Luminosity::Luminosity(bool readBunchInfo, bool readDataInfo)
{
  if(readBunchInfo)
  {
    readBunchlist();   // from totem
    readRuns();        // from cmswbm
    readBxbyBxLumi();  // from lumi group
    readEfficiency();  // based on MC
  }

  if(readDataInfo)
    readTakenEvents(); // from processed data

  readExceptions();
}

/*****************************************************************************/
Luminosity::~Luminosity()
{
}

/*****************************************************************************/
void Luminosity::readBunchlist()
{
  cerr << Helper::col(2) << " reading injection schemes" << Helper::col();
  ifstream fileI("../lumi/injectionSchemes.txt");

  string bunchSpacing, nBunches, modifier;

  while(fileI >> bunchSpacing >> nBunches >> modifier)
  {
    string injSch;
    if(modifier != "--") injSch = bunchSpacing +"_"+ nBunches +"_"+ modifier;
                    else injSch = bunchSpacing +"_"+ nBunches;

    // read filled and selected bxs
    ifstream fileB("../lumi/" + injSch +".txt");

    int bx;
    while(fileB >> bx)
    {
      bunchCrossing[injSch][bx]++;
    }
    fileB.close();

    int nfilled=0, nselect=0;
    for(auto & b : bunchCrossing[injSch])
    {
      nfilled++;

      if(b.second == selected) nselect++;
    }

    nBunches_diag[injSch] = nselect;
    nBunches_para[injSch] = nfilled;

    cerr << ".";
  }

  fileI.close();
 
  cerr << " [done]" << endl;
}

/*****************************************************************************/
string Luminosity::getInjSch(int run, int ls)
{
  string injSch;

  for(auto & runScheme : runSchemes[run])
  if(ls >= runScheme.fromLs)
    injSch = runScheme.injSch;

  return injSch;
}

//
int Luminosity::getBunches_diag(int run, int ls)
{
  string injSch = getInjSch(run, ls);

  return nBunches_diag[injSch];
}

//
int Luminosity::getBunches_para(int run, int ls)
{
  string injSch = getInjSch(run, ls); 

  return nBunches_para[injSch];
}

/*****************************************************************************/
void Luminosity::readBxbyBxLumi()
{
  cerr << Helper::col(2) << " reading bx-by-bx lumi" << Helper::col();
  igzstream fileD("../lumi/lumi_DCSONLY.csv.gz");

  int n = 0;
  string line;
  while(getline(fileD, line))
  {
    const vector<string> arg = Helper::parse(line,",");

    if(arg[0].front() != '#')
    {
      const vector<string> a = Helper::parse(arg[0],":");

      const int run  = stoi(a[0]);
      const int fill = stoi(a[1]);

      fillInfo[run] = fill;

      const int ls   = stoi(Helper::parse(arg[1],":")[0]);

      const string time   = arg[2];
      const string source = arg[8];

      //
      string injSch = getInjSch(run, ls);

      //
      LumiInfo info;

      info.time   = time;
      info.source = source;

      info.Lpara = 0;
      info.Ldiag = 0;
 
      // 
      const vector<string> list = Helper::parse(arg[9].substr(1, arg[9].size()-2)," ");

      const int nbx = list.size() / 3;

      for(int i = 0; i < nbx; i++)
      {
        const int bx = stoi(list[3*i+0]);

        if(bunchCrossing[injSch].count(bx) > 0)
        {

//        const double delLumi = stof(list[3*i+1]) / 1e-6;
          const double recLumi = stof(list[3*i+2]) / 1e-6;

          int nTimes = bunchCrossing[injSch].at(bx);

          if(nTimes >= filled  ) info.Lpara += recLumi;
          if(nTimes == selected) info.Ldiag += recLumi;
        }
      }

      lumiInfo[run][ls] = info; // store
    }

    if(++n % 250 == 0) cerr << ".";
  }

  fileD.close();

  cerr << " [done]" << endl;
}

/*****************************************************************************/
void Luminosity::readRuns()
{
  cerr << Helper::col(2) << " reading runs" << Helper::col();
  ifstream fileR("../lumi/runs.dat");

  int run; // fill, nBunches;
  string scheme;

  string line;
  while(getline(fileR, line))
  {
    const vector<string> arg = Helper::parse(line," ");

//  if(arg[0] != "s") fill     = stoi(arg[0]);
    if(arg[1] != "s") run      = stoi(arg[1]);
//  if(arg[2] != "s") nBunches = stoi(arg[2]);
    if(arg[3] != "s") scheme   =      arg[3] ;

    for(size_t i = 4; i < arg.size(); i++)
    {
      const vector<string> list = Helper::parse(arg[i].substr(1,arg[i].size()-2),",");

      RunScheme info;

      info.fromLs = stoi(list[0]);

      if(list[1] != "") info.injSch = scheme + list[1];
                   else info.injSch = scheme;

      runSchemes[run].push_back(info);
    }
    cerr << ".";
  }

  fileR.close();

  cerr << " [done]" << endl;
}

/*****************************************************************************/
void Luminosity::readTakenEvents()
{
  cerr << Helper::col(2) << " reading taken events per ls" << Helper::col();

  ifstream fileT("../out/lumi/takenEvents.dat"); 

  int run,ls;
  double tb,bt, tt,bb;

  int n = 0;
  while(fileT >> run >> ls >> tb >> bt >> tt >> bb)
  {
    takenEvents[run][ls] = {tb,bt,tt,bb};

    if(++n % 500 == 0) cerr << ".";
  }

  fileT.close();

  cerr << " [done]" << endl;
}

/*****************************************************************************/
void Luminosity::readEfficiency()
{
  cerr << Helper::col(2) << " reading average roman pot efficiency"
       << Helper::col();

  ifstream fileE("../pars/mc_efficiency.par");

  string line;
  int topo = 0;
  while(getline(fileE, line))
  {
    const vector<string> arg = Helper::parse(line," ");

    efficiency[topo] = stof(arg[ 2])*1e-2;
 
    topo++;

    cerr << ".";
  } 

  fileE.close();

  cerr << " [done]" << endl;
}

/*****************************************************************************/
void Luminosity::readExceptions()
{
  cerr << Helper::col(2) << " reading lumi exceptions" << Helper::col();

  ifstream fileE("../lumi/lumi_exceptions.par");

  string  run,lsrange, conf, stat;
  string orun,        oconf,ostat;

  int i = 0;

  while(fileE >> run >> lsrange >> conf >> stat)
  {
    if(run  == "s") run  = orun;
    if(conf == "s") conf = oconf;
    if(stat == "s") stat = ostat;

    vector<string> ls = Helper::parse(lsrange,"-");

    vector<int> topos;
    if(conf == "all" ) topos = {TB,BT,TT,BB}; 
    if(conf == "diag") topos = {TB,BT}; 
    if(conf == "para") topos = {TT,BB}; 
    if(conf == "1B"  ) topos = {BT,BB};
    if(conf == "2T"  ) topos = {BT,TT}; 

    int status = 0;
    if(stat == "off"    ) status = 0;
    if(stat == "partial") status = 1;
    if(stat == "low"    ) status = 2;

    for(auto & topo : topos)
    {
      if(ls.size() == 1)
      {
          int j = stoi(ls[0]);
          exception[stoi(run)][j][topo] = status;
      }
      else
      {
        for(int j = stoi(ls[0]); j <= stoi(ls[1]); j++)
          exception[stoi(run)][j][topo] = status;
      }
    }

    orun = run; oconf = conf; ostat = stat;

    if(++i % 10 == 0) cerr << ".";
  }

  fileE.close();

  cerr << " [done]" << endl;
}

/*****************************************************************************/
bool Luminosity::getStatus(int run, int ls, int topo, int & status)
{
  if(exception.count(run) > 0)
  if(exception[run].count(ls) > 0)
  if(exception[run][ls].count(topo) > 0)
  {
    status = exception[run][ls][topo];
    return true;
  }
 
  status = -1;
  return false;
}

/*****************************************************************************/
bool Luminosity::noProblem(int run, int ls)
{
  bool ok = true;

  if(exception.count(run) > 0)
  if(exception[run].count(ls) > 0)
    for(int topo = 0; topo < nTopos; topo++)
      if(exception[run][ls].count(topo) > 0) ok = false;

  return ok;
}

/*****************************************************************************/
// average roman pot efficiency
double Luminosity::getEfficiency(int topo) 
{
  return efficiency[topo];
}

/*****************************************************************************/
// special efficiency
double Luminosity::specEfficiency(int run, int ls, int topo)
{
  double eff = 1;

  //
  if(run == 319159 && ls <  250)
                                 eff *= 437./ 732;

  if(run == 319222 && ls >= 233)
    if(topo == TB || topo == BT) eff *= 244./ 487;

  if(run == 319300 && ls <  210)
    if(topo == TB || topo == BT) eff *= 483./1450;

  return eff;
}

/*****************************************************************************/
// recording efficiency
double Luminosity::recEfficiency(int status)
{
  //
  const double lowEfficiency = 3./4; // PARAMETER

  return (status ==-1 ? 1 : // normal
         (status == 0 ? 0 : // off
         (status == 1 ? 1 : // partial
          lowEfficiency))); // low
}

/*****************************************************************************/
// average collisions in a bunch crossing
double Luminosity::getMu(double Lpara, int nBunches, double sig_vis)
{
  const double nOrbits = pow(2,18); // 2^18, means 23.3 sec

  return Lpara/(nBunches*nOrbits) * sig_vis;
}

/*****************************************************************************/
// predicted
double Luminosity::getEv(double L, double mu, double sig)
{
  return L*sig * exp(-mu);
}

/*****************************************************************************/
void Luminosity::calcWeights()
{
  //
  cerr << " processing events";

  vector<double>           recLumi_all(nTopos,0); // recLumi_all[topo]
  map<int,map<int,double>> recLumi_run;           // recLumi_run[run][topo]

  vector<double>           effLumi_all(nTopos,0); // effLumi_all[topo]
  map<int,map<int,double>> effLumi_run;           // effLumi_run[run][topo]

  ofstream fileL("../out/lumi/lumi.dat");

  //
  int i = 0;
  for(auto & a : lumiInfo)
  {
    const int & run = a.first;

    for(auto & b : a.second)
    {
      const int & ls = b.first;

      const LumiInfo & info = b.second;

      // diag
      double Ldiag = info.Ldiag;
      double Lpara = info.Lpara;

      const int nBunches_diag = getBunches_diag(run, ls);
      const int nBunches_para = getBunches_para(run, ls);

      const double mu_diag = getMu(Ldiag, nBunches_diag, sig_vis);
      const double mu_para = getMu(Lpara, nBunches_para, sig_vis);

      // fix run 319260
      if(run == 319260 && ls >= 30)
      {
        Ldiag = (-0.0303*ls + 520.9435) / 1e-6;
        Lpara = (-0.0429*ls + 676.9360) / 1e-6;
      }

      //
      fileL << " " << run << " " << ls;

      for(int topo = 0; topo < nTopos; topo++)
      {
        int status;
        getStatus(run,ls, topo,status);

        // recorded
        double Lrec = (topo == TB || topo == BT ?   Ldiag :   Lpara) *
                       specEfficiency(run,ls,topo) *
                        recEfficiency(status);

        // pileup
        double mu   = (topo == TB || topo == BT ? mu_diag : mu_para);

        // effective
        double Leff = Lrec * exp(-mu);

        fileL << " " << Leff;

        recLumi_all[topo]      += Lrec;
        recLumi_run[run][topo] += Lrec;

        effLumi_all[topo]      += Leff;
        effLumi_run[run][topo] += Leff;
      }

      fileL << " " << info.time << endl;

      //
      if(++i % 1000 == 0) cerr << ".";
    }
  }

  fileL.close();

  cerr << " [done]" << endl;

  //
  ofstream fileE("../out/lumi/int_lumi.dat");
  cerr << " intLumi:" << endl;
  for(int topo = 0; topo < nTopos; topo++)
  {
    cerr  << " " << topos[topo]
          << " " << int(recLumi_all[topo]*1e-12 * 1e+3)/1e+3
          << " " << int(effLumi_all[topo]*1e-12 * 1e+3)/1e+3
          << " [pb-1]"
          << endl;

    fileE << " " << topos[topo] << " & " << recLumi_all[topo]
                                << " & " << effLumi_all[topo] << endl;
  }
  cerr << " [pb^-1]" << endl;
  fileE.close();

  //
  ofstream fileR("../out/lumi/run_lumi.dat");
  for(auto & m : recLumi_run)
  {
    const int & run = m.first;
    const map<int,double> & val = m.second;

    fileR << " " << run;

    for(int topo = 0; topo < nTopos; topo+=2)
    {
      double v = (val.at(topo) +
                  val.at(topo+1))/2.;

      fileR << " " << int(v*1e-12 * 1000)/1000.;
    }

    fileR << endl; // [pb^{-1}]
  }
  fileR.close();
}

/*****************************************************************************/
void Luminosity::postProcess()
{
  //
  cerr << " processing events";

  ofstream fileO("../out/lumi/lumi.out");

  int i = 0;
  for(auto & a : takenEvents)
  {
    const int & run = a.first;

    for(auto & b : a.second)
    {
      const int & ls = b.first;

      const vector<double> & measEvents = b.second;
      const LumiInfo & info = lumiInfo[run][ls];

      vector<double> predEvents(nTopos);

      //
      double Ldiag = info.Ldiag;
      double Lpara = info.Lpara;

      const int nBunches_diag = getBunches_diag(run, ls);
      const int nBunches_para = getBunches_para(run, ls);

      const double mu_diag = getMu(Ldiag, nBunches_diag, sig_vis);
      const double mu_para = getMu(Lpara, nBunches_para, sig_vis);

      // fix run 319260
      if(run == 319260 && ls >= 30)
      {
        Ldiag = (-0.0303*ls + 520.9435) / 1e-6;
        Lpara = (-0.0429*ls + 676.9360) / 1e-6;
      }

      
      { // diag
        const double nEvents = getEv(Ldiag, mu_diag, sig_diag);

        predEvents[TB] = nEvents * getEfficiency(TB);
        predEvents[BT] = nEvents * getEfficiency(BT);
      }

      {
        // para
        const double nEvents = getEv(Lpara, mu_para, sig_para);

        predEvents[TT] = nEvents * getEfficiency(TT);
        predEvents[BB] = nEvents * getEfficiency(BB);
      }

      //
      vector<int> status(nTopos);
      for(int topo = 0; topo < nTopos; topo++)
      {
        getStatus(run,ls,topo, status[topo]);
        predEvents[topo] *= specEfficiency(run,ls,topo);
      }

      //
      fileO << " " << predEvents[TB] << " " << measEvents[TB]
            << " " << predEvents[BT] << " " << measEvents[BT]
            << " " << predEvents[TT] << " " << measEvents[TT]
            << " " << predEvents[BB] << " " << measEvents[BB] << " "
            << " " << fillInfo[run] //  9
            << " " << run           // 10
            << " " << ls            // 11
            << " " << info.time     // 12
            << " " << status[TB] // 14
            << " " << status[BT] // 15
            << " " << status[TT] // 16
            << " " << status[BB] // 17
            << " " << mu_diag // 18
            << " " << mu_para // 19
            << endl;

      //
      if(++i % 1000 == 0) cerr << ".";
    }

    fileO << endl;
  }

  fileO.close();

  cerr << " [done]" << endl;
}

