#include "../interface/HistoVal.h"

#include "../interface/Random.h"
#include "../interface/gzstream.h"

#include <cmath>
#include <algorithm>

#define sqr(x) ((x)*(x))

using namespace std;

/*****************************************************************************/
HistoVal::HistoVal()
{
}

/*****************************************************************************/
HistoVal::~HistoVal()
{
  if(toWrite)
    write();
}

/*****************************************************************************/
float HistoVal::getVolume()
{
  float vol = 1;

  for(auto & axis : axes)
    vol *= axis.binwidth();

  return vol;
}

/*****************************************************************************/
void HistoVal::init(float xmin, float xmax, int xbins,
                    float ymin, float ymax, int ybins,
                    float zmin, float zmax, int zbins,
                    float vmin, float vmax, int vbins,
                    string name_, bool norm_)
{
  type = 4; name = name_; norm = norm_;
  axes.push_back({xmin,xmax,xbins});
  axes.push_back({ymin,ymax,ybins});
  axes.push_back({zmin,zmax,zbins});
  axes.push_back({vmin,vmax,vbins});

  a4.clear(); a4.resize(xbins);
  for(int ix = 0; ix < xbins; ix++)
  {
    a4[ix].resize(ybins);

    for(int iy = 0; iy < ybins; iy++)
    {
      a4[ix][iy].resize(zbins);

      for(int iz = 0; iz < zbins; iz++)
        a4[ix][iy][iz].resize(vbins, 0);
    }
  }

  binvolume = getVolume();
}

/*****************************************************************************/
bool HistoVal::index(const vector<float> & x, vector<int> & ix)
{
  bool ok = true;

  for(size_t i = 0; i < x.size(); i++)
  {
    const auto & axis = axes[i];

    if(!(x[i] >= axis.min && x[i] < axis.max))
    {
      ok = false; break;
    }
  }

  if(ok)
    for(size_t i = 0; i < axes.size(); i++)
    {
      const auto & axis = axes[i];

      int j = int( (x[i] - axis.min) / (axis.max - axis.min) * axis.bins);

      if(j < 0)          j = 0;
      if(j >= axis.bins) j = axis.bins - 1;

      ix.push_back(j);
    }

  return ok;
}

/*****************************************************************************/
void HistoVal::add(float & a, const float & v)
{
  mtx.lock();   // MUTEX

  a += v;   // num

  mtx.unlock(); // MUTEX
}

/*****************************************************************************/
float & HistoVal::get(const vector<int> & ix)
{
  switch(type)
  {
    case 4 :  return               a4[ix[0]][ix[1]][ix[2]].at(ix[3]); break;
    default : exit(1);
  }

  cerr << " Histo::get (" << name << ") problem" << endl; exit(1);
}

/*****************************************************************************/
float HistoVal::val(const vector<float> & x)
{
  vector<int> ix;

  if(index(x, ix))
    return get(ix);
  else
    return 0;
}

/*****************************************************************************/
void HistoVal::fillw(const vector<float> & x, float w)
{
  vector<int> ix;

  if(int(x.size()) != type)
  {
    cerr << " Histo::fillw (" << name << ") problem : x.size() != type" << endl;
    exit(1);
  }

  if(index(x, ix))
    add(get(ix), w);
}

//
void HistoVal::fill(const vector<float> & x)
{
  fillw(x, 1.);
}

/*****************************************************************************/
void HistoVal::set(const vector<float> & x, float v)
{
  vector<int> ix;

  if(index(x, ix))
    get(ix) = v;
}

/*****************************************************************************/
void HistoVal::div(float & a, const float & b)
{
        float & vala = a;
  const float & valb = b;

  float val = (valb != 0 ? vala/valb : 0); 

  vala = val;
}

//
void HistoVal::div(const HistoVal & q)
{
  if(type == 4)
    for(int ix = 0; ix < axes[0].bins; ix++)
    for(int iy = 0; iy < axes[1].bins; iy++)
    for(int iz = 0; iz < axes[2].bins; iz++)
    for(int iv = 0; iv < axes[3].bins; iv++)
      div(a4[ix][iy][iz][iv], q.a4[ix][iy][iz][iv]);
}

/*****************************************************************************/
void HistoVal::write(ofstream & file, vector<int> & ix, vector<float> & x, int j)
{
  if(j < type)
  {
    const auto & axis = axes[j];

    if(norm && j == type - 1)
    {
      sum = 0;
      for(ix[j] = 0; ix[j] < axis.bins; ix[j]++)
        sum += get(ix);

      sum *= (axis.max - axis.min) / axis.bins;
    }

    for(ix[j] = 0; ix[j] < axis.bins; ix[j]++)
    {
      x[j] = axis.min + (ix[j]+0.5)/axis.bins * (axis.max - axis.min);
      write(file, ix,x, j+1);
    }

    if(j == type - 2) 
      file << endl;

    if(j == type - 1) 
      file << endl;
  }
  else
  {
    for(int k = 0; k < type; k++)
      file << " " << x[k];

    float elem = get(ix);

    if(norm && sum != 0)
      file << " " <<      elem   / sum
           << endl;   
    else
    {
      file << " " <<      elem        / (byVol ? binvolume : 1)
           << endl;
    }
  }
}

//
void HistoVal::write()
{
  if(name == "") return;

  ofstream file(name);

  vector<int>   ix(type,0);
  vector<float> x(type,0);

  sum = 0;
  write(file, ix,x, 0);

  file << endl << endl;

  file.close();
}

/*****************************************************************************/
void HistoVal::read(igzstream & file, vector<int> & ix, int j)
{
  if(j < type)
  {
    const auto & axis = axes[j];

    for(ix[j] = 0; ix[j] < axis.bins; ix[j]++)
      read(file, ix, j+1);
  }
  else
  {
    float f;
    for(int k = 0; k < type; k++)
      file >> f;

    float val;

    file >> val;

    float & elem = get(ix);

    elem = val;
  }
}

//
void HistoVal::read(ifstream & file, vector<int> & ix, int j)
{
  if(j < type)
  {
    const auto & axis = axes[j];

    for(ix[j] = 0; ix[j] < axis.bins; ix[j]++)
      read(file, ix, j+1);
  }
  else
  {
    float f;
    for(int k = 0; k < type; k++)
      file >> f;

    float val;

    file >> val;

    float & elem = get(ix);

    elem = val;
  }
}

//
void HistoVal::read(bool fromZip)
{
  if(fromZip)
  {
    igzstream file(name.c_str());
    vector<int> ix(type,0);
    read(file, ix, 0);
    file.close();
  }
  else
  {
    ifstream file(name);
    vector<int> ix(type,0);
    read(file, ix, 0);
    file.close();
  }

  toWrite = false;
}

/*****************************************************************************/
void HistoVal::normalize()
{
  if(type == 3)
  {
    vector<int> ix(3);

    for(ix[0] = 0; ix[0] < axes[0].bins; ix[0]++)
    for(ix[1] = 0; ix[1] < axes[1].bins; ix[1]++)
    {
      float sum = 0;
      for(ix[2] = 0; ix[2] < axes[2].bins; ix[2]++)
        sum += get(ix);

      float cum = 0;
      for(ix[2] = 0; ix[2] < axes[2].bins; ix[2]++)
      {
        cum += get(ix);

        if(sum > 0) get(ix) = cum / sum;
               else get(ix) = -1;
      }
    }
  }
  else
  { cerr << " Histo::normalize (" << name << ") problem" << endl; exit(1); }
}

/*****************************************************************************/
void HistoVal::normalizeTrue()
{
  if(type == 2)
  {
    vector<int> ix(2);

    for(ix[0] = 0; ix[0] < axes[0].bins; ix[0]++)
    {
      float sum = 0;
      for(ix[1] = 0; ix[1] < axes[1].bins; ix[1]++)
        sum += get(ix);

      for(ix[1] = 0; ix[1] < axes[1].bins; ix[1]++)
      {
        if(sum > 0)
        {
          get(ix)  /= sum;
        }
        else
          get(ix) = -1;
      }
    }
  }

  if(type == 3)
  {
    vector<int> ix(3);

    for(ix[0] = 0; ix[0] < axes[0].bins; ix[0]++)
    for(ix[1] = 0; ix[1] < axes[1].bins; ix[1]++)
    {
      float sum = 0;
      for(ix[2] = 0; ix[2] < axes[2].bins; ix[2]++)
        sum += get(ix);

      for(ix[2] = 0; ix[2] < axes[2].bins; ix[2]++)
      {
        if(sum > 0)
        {
          get(ix) /= sum;
        }
        else
          get(ix) = -1;
      }
    }
  }
}

