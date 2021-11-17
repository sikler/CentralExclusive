#ifndef _RpSimu_h_
#define _RpSimu_h_

#include <map>
#include <vector>
#include <mutex>

//
struct RpTrack;
struct RpDet;

class RpFit;
class RpEffic;
class Histo;
class Random;
class ogzstream;

//
class RpSimu
{
 public:
  RpSimu(int nCpus, int nEvents);
  ~RpSimu();

  void generateData(int nCmTra, int run);

 private:
   bool isSingle(double x);
   double generateHit();

   bool generateRpTrack(int run, RpTrack & track);
   void printRpTrack(ogzstream & file, RpTrack & track);
   void generateStream(int run, const std::vector<std::string> & dataStreams);

   RpFit   * theRpFit;
   RpEffic * theRpEffic;
   Random  * theRandom;

   std::map<RpDet, Histo> his_effic, den_effic;

   std::mutex mtx_rnd, mtx_his, mtx_count;

   const int nCpus, nEvents;
   int iEvents;

   bool toFile;
};

#endif
