#ifndef _HistoVal_h_
#define _HistoVal_h_

#include <fstream>
#include <string>
#include <vector>
#include <mutex>

class igzstream;
class Random;

#ifndef _Axis_
#define _Axis_
struct Axis {
  float min,max;
  int bins;

  float center(int ix) const
  { return min + (max - min)/bins * (ix+0.5); }

  float binwidth() const
  { return (max - min)/bins; }
};
#endif

class HistoVal
{
 public:
  HistoVal();
  ~HistoVal(); 

  void init(float xmin, float xmax, int xbins,
            float ymin, float ymax, int ybins,
            float zmin, float zmax, int zbins,
            float vmin, float vmax, int vbins,
            std::string name, bool norm = false);

  float val(const std::vector<float> & x);

  void fill (const std::vector<float> & x);
  void fillw(const std::vector<float> & x, float w);

  void set(const std::vector<float> & x, float v);

  void div(const HistoVal & q);

  void write();
  void read(bool fromZip=false);
 
  void normalize();
  void normalizeTrue();

  //
  float & get(const std::vector<int> & ix);
  std::vector<Axis> axes;

  bool byVol = false; // divide by volume (for physics)

  bool toWrite = true; // FIXME

 private:
  HistoVal(const HistoVal &);
  HistoVal(HistoVal &);

  float getVolume();

  void add(float & a, const float & v);

  bool index(const std::vector<float> & x, std::vector<int> & ix);

  void div(float & a, const float & b);

  void write(std::ofstream & file,
             std::vector<int> & ix, std::vector<float> & x, int j);

  void read (igzstream & file,
             std::vector<int> & ix, int j);
  void read (std::ifstream & file,
             std::vector<int> & ix, int j);

  // (val,sig2)
  std::vector<std::vector<std::vector<std::vector<float>>>> a4;

  float binvolume;

  int type;

  std::string name;
  bool norm;

  float sum;

  std::mutex mtx;
};

#endif

