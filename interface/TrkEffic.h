#ifndef _TrkEffic_h_
#define _TrkEffic_h_

#include <vector>
#include <map>

#include "../interface/Histo.h"

#include "../interface/Parameters.h"

//
struct Part;
struct SimTrackInfo;

struct CmTrack;
struct Vector3;

class Resolutions;

//
class TrkEffic
{
 public:
  TrkEffic(int flag);
  ~TrkEffic();

  void processSimulation(const std::string & partName);

  double getEfficiency(const std::vector<CmTrack> & tracks,    int type);
  double getEfficiency(const Vector3 & p3, const Vector3 & p4, int type);

 private:
  int getIndex(double low, double hig, int nbins, double x);
  bool getIndex(const Part & p, std::vector<short int> & index);

  void readSimTable(const std::string & partName,
                    std::map<std::vector<short int>,SimTrackInfo> & simInfo);
  void readSimTables();

  int getLayerCountBpix(int layerCodes);

  bool isCompatible(const Part & gen, const Part & rec);
 
  double getEfficiency(Part p[2], int type);

  //
  std::map<std::vector<short int>,SimTrackInfo> simInfo[nParts][2]; // charge

  Resolutions * theResolutions;
};

#endif
