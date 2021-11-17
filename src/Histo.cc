#include "../interface/Histo.h"

#include "../interface/Random.h"
#include "../interface/gzstream.h"

#include <cmath>
#include <algorithm>

#define sqr(x) ((x)*(x))

using namespace std;

/*****************************************************************************/
Histo::Histo()
{
}

/*****************************************************************************/
Histo::~Histo()
{
  if(toWrite)
    write();
}

/*****************************************************************************/
float Histo::getVolume()
{
  float vol = 1;

  for(auto & axis : axes)
    vol *= axis.binwidth();

  return vol;
}

/*****************************************************************************/
void Histo::init(float xmin, float xmax, int xbins,
                 string name_, bool norm_)
{
  type = 1; name = name_; norm = norm_;
  axes.push_back({xmin,xmax,xbins});

  a1.clear(); a1.resize(xbins, pair<float,float>(0,0));

  binvolume = getVolume();
}

//
void Histo::init(float xmin, float xmax, int xbins,
                 float ymin, float ymax, int ybins,
                 string name_, bool norm_)
{
  type = 2; name = name_; norm = norm_;
  axes.push_back({xmin,xmax,xbins});
  axes.push_back({ymin,ymax,ybins});

  a2.clear(); a2.resize(xbins);
  for(int ix = 0; ix < xbins; ix++)
    a2[ix].resize(ybins, pair<float,float>(0,0));

  binvolume = getVolume();
}

//
void Histo::init(float xmin, float xmax, int xbins,
                 float ymin, float ymax, int ybins,
                 float zmin, float zmax, int zbins,
                 string name_, bool norm_)
{
  type = 3; name = name_; norm = norm_;
  axes.push_back({xmin,xmax,xbins});
  axes.push_back({ymin,ymax,ybins});
  axes.push_back({zmin,zmax,zbins});

  a3.clear(); a3.resize(xbins);
  for(int ix = 0; ix < xbins; ix++)
  {
    a3[ix].resize(ybins);

    for(int iy = 0; iy < ybins; iy++)
      a3[ix][iy].resize(zbins, pair<float,float>(0,0));
  }

  binvolume = getVolume();
}

//
void Histo::init(float xmin, float xmax, int xbins,
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
        a4[ix][iy][iz].resize(vbins, pair<float,float>(0,0));
    }
  }

  binvolume = getVolume();
}

/*****************************************************************************/
bool Histo::index(const vector<float> & x, vector<int> & ix)
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
void Histo::add(pair<float,float> & a, const float & v)
{
  mtx.lock();   // MUTEX

  a.first  += v;   // num
  a.second += v*v; // sig^2

  mtx.unlock(); // MUTEX
}

/*****************************************************************************/
pair<float,float> & Histo::get(const vector<int> & ix)
{
  switch(type)
  {
    case 1 : return                                    a1.at(ix[0]); break;
    case 2 : return                             a2[ix[0]].at(ix[1]); break;
    case 3 : return                      a3[ix[0]][ix[1]].at(ix[2]); break;
    case 4 : return               a4[ix[0]][ix[1]][ix[2]].at(ix[3]); break;
    case 5 : exit(1); break;
    case 6 : exit(1); break;
  }

  cerr << " Histo::get (" << name << ") problem" << endl; exit(1);
}

/*****************************************************************************/
pair<float,float> Histo::getv(const vector<float> & x)
{
  vector<int> ix;

  if(index(x, ix))
    return get(ix);
  else
    return pair<float,float>(0,0);
}

/*****************************************************************************/
float Histo::val(const vector<float> & x)
{
  vector<int> ix;

  if(index(x, ix))
    return get(ix).first;
  else
    return 0;
}

/*****************************************************************************/
void Histo::fillw(const vector<float> & x, float w)
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
void Histo::fill(const vector<float> & x)
{
  fillw(x, 1.);
}

/*****************************************************************************/
void Histo::set(const vector<float> & x, float v)
{
  vector<int> ix;

  if(index(x, ix))
    get(ix) = pair<float,float>(v, sqr(v));
}

/*****************************************************************************/
void Histo::div(pair<float,float> & a, const pair<float,float> & b)
{
        float & vala = a.first;
  const float & valb = b.first;

        float & si2a = a.second;
  const float & si2b = b.second;

  float val = (valb != 0 ? vala/valb : 0); 
  float si2 = sqr(val) * (si2a/sqr(vala) + si2b/sqr(valb));

  vala = val;
  si2a = si2;
}

//
void Histo::div(const Histo & q)
{
  if(type == 1)
    for(int ix = 0; ix < axes[0].bins; ix++)
      div(a1[ix], q.a1[ix]);

  if(type == 2)
    for(int ix = 0; ix < axes[0].bins; ix++)
    for(int iy = 0; iy < axes[1].bins; iy++)
      div(a2[ix][iy], q.a2[ix][iy]);

  if(type == 3)
    for(int ix = 0; ix < axes[0].bins; ix++)
    for(int iy = 0; iy < axes[1].bins; iy++)
    for(int iz = 0; iz < axes[2].bins; iz++)
      div(a3[ix][iy][iz], q.a3[ix][iy][iz]);

  if(type == 4)
    for(int ix = 0; ix < axes[0].bins; ix++)
    for(int iy = 0; iy < axes[1].bins; iy++)
    for(int iz = 0; iz < axes[2].bins; iz++)
    for(int iv = 0; iv < axes[3].bins; iv++)
      div(a4[ix][iy][iz][iv], q.a4[ix][iy][iz][iv]);

  if(type == 5) exit(1);

  if(type == 6) exit(1);
}

/*****************************************************************************/
void Histo::write(ofstream & file, vector<int> & ix, vector<float> & x, int j)
{
  if(j < type)
  {
    const auto & axis = axes[j];

    if(norm && j == type - 1)
    {
      sum = 0;
      for(ix[j] = 0; ix[j] < axis.bins; ix[j]++)
        sum += get(ix).first;

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

    pair<float,float> elem = get(ix);

    if(norm && sum != 0)
      file << " " <<      elem.first   / sum
           << " " << sqrt(elem.second > 0 ? elem.second : 1) / sum
           << endl;   
    else
    {
      file << " " <<      elem.first        / (byVol ? binvolume : 1)
           << " " << sqrt(elem.second > 0 ?
                          elem.second : 1)  / (byVol ? binvolume : 1)
           << endl;
    }
  }
}

//
void Histo::write()
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
void Histo::read(igzstream & file, vector<int> & ix, int j)
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

    float val,sig = 0;

    file >> val >> sig;

    pair<float,float> & elem = get(ix);

    elem = pair<float,float>(val,sqr(sig));
  }
}

//
void Histo::read(ifstream & file, vector<int> & ix, int j)
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

    float val,sig = 0;

    file >> val >> sig;

    pair<float,float> & elem = get(ix);

    elem = pair<float,float>(val,sqr(sig));
  }
}

//
void Histo::read(bool fromZip)
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
void Histo::normalize()
{
  if(type == 3)
  {
    vector<int> ix(3);

    for(ix[0] = 0; ix[0] < axes[0].bins; ix[0]++)
    for(ix[1] = 0; ix[1] < axes[1].bins; ix[1]++)
    {
      float sum = 0;
      for(ix[2] = 0; ix[2] < axes[2].bins; ix[2]++)
        sum += get(ix).first;

      float cum = 0;
      for(ix[2] = 0; ix[2] < axes[2].bins; ix[2]++)
      {
        cum += get(ix).first;

        if(sum > 0) get(ix).first = cum / sum;
               else get(ix).first = -1;
      }
    }
  }
  else
  { cerr << " Histo::normalize (" << name << ") problem" << endl; exit(1); }
}

/*****************************************************************************/
void Histo::normalizeTrue()
{
  if(type == 2)
  {
    vector<int> ix(2);

    for(ix[0] = 0; ix[0] < axes[0].bins; ix[0]++)
    {
      float sum = 0;
      for(ix[1] = 0; ix[1] < axes[1].bins; ix[1]++)
        sum += get(ix).first;

      for(ix[1] = 0; ix[1] < axes[1].bins; ix[1]++)
      {
        if(sum > 0)
        {
          get(ix).first  /= sum;
          get(ix).second /= sum;
        }
        else
          get(ix).first = -1;
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
        sum += get(ix).first;

      for(ix[2] = 0; ix[2] < axes[2].bins; ix[2]++)
      {
        if(sum > 0)
        {
          get(ix).first  /= sum;
          get(ix).second /= sum;
        }
        else
          get(ix).first = -1;
      }
    }
  }
}

/*****************************************************************************/
float Histo::sample(const vector<float> & x, Random * theRandom)
{
  vector<int> ix;
  if(index(x, ix))
  {
    const vector<pair<float,float>> & v = a3[ix[0]][ix[1]];
    
    float y = theRandom->getFlat();   

    // binary search
    vector<pair<float,float>>::const_iterator up =
      upper_bound(v.begin(), v.end(), pair<float,float>(y,0));

    if(up == v.begin())
    { cerr << " Histo::sample (" << name << ") problem" << endl; exit(1); }

    const float dz = (axes[2].max - axes[2].min)/axes[2].bins;

    float x_up = axes[2].min + (up - v.begin()) * dz;

    float y_up =     up->first;
    float y_lo = (up-1)->first;

    float r = x_up + (y_up - y)/(y_up - y_lo) * dz;

    return r;
  }
  else
    return 0;
}

