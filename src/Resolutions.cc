#include "../interface/Resolutions.h"

#include "../interface/Structures.h"
#include "../interface/Random.h"

#include <iostream>
#include <string>

#define sqr(x) ((x)*(x))

using namespace std;

// flag = 1 : pt resolution // simu
// flag = 2 : mass response // doPhys

/*****************************************************************************/
Resolutions::Resolutions(int flag, const string & partName) : flag(flag) // for simu, pt resolution
{
  ptResol.init(-maxEta,maxEta,nEta,
                     0,maxPt ,nPt ,
               -0.5,0.5,200,
               "../out/track/ptResol_"+partName+".his"); // ,true); don't norm
}

/*****************************************************************************/
Resolutions::Resolutions(int flag) : flag(flag) // for data, mass resolution
{
  theRandom = new Random();

  for(int type = 0; type < nParts; type++)
  for(int q = 0; q < 2; q++)
  {
    { Histo & h = ptResol_mean[type][q];
      h.init(-maxEta,maxEta,nEta, 0,maxPt,nPt, "");
      h.toWrite = false; }

    { Histo & h = ptResol_sigm[type][q];
      h.init(-maxEta,maxEta,nEta, 0,maxPt,nPt, "");
      h.toWrite = false; }
  }

  //
  cerr << " reading pt resolution..";
  for(int type = 0; type < nParts; type++)
  for(int q = 0; q < 2; q++)
  {
    ifstream file("../out/track/ptResol_"+parts[type]+(q==0?"p":"m")+".dat");

    while(!file.eof())
    {
      double eta,pt, mean,sigm;
      file >> eta >> pt >> mean >> sigm;

      if(!file.eof())
      {
        ptResol_mean[type][q].set({eta,pt}, mean);
        ptResol_sigm[type][q].set({eta,pt}, sigm);

        if(pt == 0.075) // copy to 25 MeV as well
        {
          ptResol_mean[type][q].set({eta,pt-0.050}, mean);
          ptResol_sigm[type][q].set({eta,pt-0.050}, sigm);
        }
      }
    }

    file.close();
  }
  cerr << " [done]" << endl;

  //
  for(int i1 = 0; i1 < nKt; i1++)
  for(int i2 = 0; i2 < nKt; i2++)
  for(int type = 0; type < nParts; type++)
  {
    Histo & h = his_massResp[i1][i2][type];
    char name[256];
    sprintf(name,"../out/track/massResol/%d_%d_%s.dat",
                               i1,i2,parts[type].c_str());

    h.init( 0,maxMass,nMass, -110e-3,110e-3, 11, name);
  }
}

/*****************************************************************************/
Resolutions::~Resolutions()
{
  if(flag == 2)
  for(int i1 = 0; i1 < nKt; i1++)
  for(int i2 = 0; i2 < nKt; i2++)
  for(int type = 0; type < nParts; type++)
  {
    Histo & h = his_massResp[i1][i2][type];
    h.normalizeTrue();
  }
}

/*****************************************************************************/
void Resolutions::collectForPt(const vector<float> & etaptdpt)
{
  ptResol.fill(etaptdpt);
}

/*****************************************************************************/
void Resolutions::collectForMass(int type, const vector<float> & dphip1tp2t,
                                 double m, double weight,
                                 const Vector4 & p3, const Vector4 & p4)
{
  // collect mass response matrix
  // write expected shift and sigma_m at m

  const double z = 1e-4; // small

  double dm_dp3t, dm_dp4t;

  {
    // increase p3 momentum vector by z (direction is fix since well known)
    Vector4 pp = {sqrt(sqr(mass[type]) + sqr((1+z)*p3.p.length())),
                 (1+z)*p3.p.x, (1+z)*p3.p.y, (1+z)*p3.p.z};

    dm_dp3t = ( (pp + p4).mass() - (p3 + p4).mass() ) / 
              (     pp.p.trans() -     p3.p.trans() );
  }

  {
    // increase p4 momentum vector by z (direction is fix since well known)
    Vector4 pp = {sqrt(sqr(mass[type]) + sqr((1+z)*p4.p.length())),
                 (1+z)*p4.p.x, (1+z)*p4.p.y, (1+z)*p4.p.z};

    dm_dp4t = ( (pp + p3).mass() - (p4 + p3).mass() ) /
                (   pp.p.trans() -     p4.p.trans() );
  }

  double shift_m =
               dm_dp3t * ptResol_mean[type][0].val({p3.p.eta(), p3.p.trans()})
             + dm_dp4t * ptResol_mean[type][1].val({p4.p.eta(), p4.p.trans()});

  double sigma_m =
      sqrt(sqr(dm_dp3t * ptResol_sigm[type][0].val({p3.p.eta(), p3.p.trans()}))
         + sqr(dm_dp4t * ptResol_sigm[type][1].val({p4.p.eta(), p4.p.trans()})));

  //
  const double & p1t  = dphip1tp2t[1];
  const double & p2t  = dphip1tp2t[2];

  const double dKt = (maxKt - minKt)/nKt;

  const int i1 = int((p1t - minKt)/dKt);
  const int i2 = int((p2t - minKt)/dKt);

  if(m < maxMass)
  if(p1t > minKt && p1t < maxKt)
  if(p2t > minKt && p2t < maxKt)
  for(int j = 0; j < 100; j++) // 100 trials FIXME
  {
    double mp = m + shift_m + theRandom->getGauss()*sigma_m;

    his_massResp[i1][i2][type].fillw({m,mp-m}, weight);
  }
}

