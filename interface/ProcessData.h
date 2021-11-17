#ifndef _ProcessData_h_
#define _ProcessData_h_

#include <vector>
#include <map>
#include <mutex>

//
struct Vector2;
struct RpDet;
struct RpTrack;
struct Event;

class RpPat;
class RpEffic;
class RpReco;
class RpVeto;

class ParticleId;
class TrkEffic;
class TrkReco;
class Selection;

class Physics;

class igzstream;

//
class ProcessData
{
 public:
  ProcessData(const std::string & data,
              int nCmTra, int nEvery, int nCpus,
              const std::vector<std::string> & flags);
  ~ProcessData();

  void readData();

 private:
  bool has(const std::string & flag);

  //
  bool readEvent(igzstream & fileIn, int nCmTra, Event & event);

  //
  bool processEvent(Event & event);
  void readStreams(const std::vector<std::string> & dataStreams);

  //
  /*const*/ std::string data;
  int nCmTra, nEvery, nCpus;
  const std::vector<std::string> flags;

  //
  std::mutex mtx_his, mtx_cerr, mtx_count;
  int nEvents;

  // [<run,ls>][topo]
  std::mutex mtx_takenMap;
  std::map<std::pair<int,int>,std::map<int,double>> takenEvents;

  //
  RpPat   * theRpPat;
  RpEffic * theRpEffic;
  RpReco  * theRpReco;
  RpVeto  * theRpVeto;

  ParticleId * theParticleId;
  TrkEffic   * theTrkEffic;
  TrkReco    * theTrkReco;
  Selection  * theSelection;

  Physics * thePhysics;
};

#endif
