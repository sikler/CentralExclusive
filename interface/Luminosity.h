#ifndef _Luminosity_h_
#define _Luminosity_h_

#include <map>
#include <vector>

//
struct LumiInfo {
  std::string time;
  std::string source;

  double Lpara;
  double Ldiag;
};

//
struct RunScheme {
  int fromLs;
  std::string injSch;
};

//
class Luminosity
{
 public:
  Luminosity();
  Luminosity(bool readBunchInfo, bool readDataInfo);
  ~Luminosity();

  void calcWeights();
  void postProcess();

  bool getStatus(int run, int ls, int topo, int & status);
  bool noProblem(int run, int ls);

 private:
  std::string getInjSch(int run, int ls);

  int   getBunches_diag(int run, int ls);
  int   getBunches_para(int run, int ls);

  void readBunchlist();
  void readBxbyBxLumi();
  void readRuns();
  void readTakenEvents();
  void readEfficiency();
  void readExceptions();

  // average roman pot efficiency
  double  getEfficiency(int topo);

  // special efficiency
  double specEfficiency(int run, int ls, int topo);

  // recording efficiency
  double  recEfficiency(int status);

  //
  double getMu(double Lpara, int nBunches, double sig);
  double getEv(double L,     double mu,    double sig);

  //
  std::map<std::string,std::map<int,int>> bunchCrossing; // [injSch][bx] = 1,2
  std::map<std::string,int> nBunches_diag; // [injSch]
  std::map<std::string,int> nBunches_para; // [injSch]

  std::map<int,std::vector<RunScheme>> runSchemes;          // [run]

  std::map<int,std::map<int,LumiInfo>> lumiInfo;            // [run][ls]
  std::map<int,std::map<int,std::vector<double>>> takenEvents; // [run][ls]

  std::map<int, double> efficiency;  // [topo]
  
  std::map<int,int> fillInfo; // [run]

  std::map<int,std::map<int,std::map<int,int>>> exception; // [run][ls][topo]

};

#endif
