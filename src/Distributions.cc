#include "../interface/Distributions.h"

#include <cmath>

#include "../interface/Vectors.h"
#include "../interface/Helper.h"

#include <iostream>

using namespace std;

/*****************************************************************************/
Distributions::Distributions()
{
  const string b = "../out/res/";

  // dphi vs dphi_tilde
  for(int type = 0; type < nParts; type++)
    his_dphis[type].init(0,maxDphi,2*nDphi,
                         0,maxDphi,2*nDphi, b + "phis_"+parts[type]+".his");

  // read mass regions
  ifstream fileReg("../pars/regions.par");
  for(int reg = 0; reg < nRegs; reg++)
    fileReg >> regType[reg] >> regName[reg] >> lowSign[reg] >> higSign[reg];
  fileReg.close();

  // dN/dphi
  for(int reg = 0; reg < nRegs; reg++)
  {
    //
    for(int topo = 0; topo < nTopos; topo++)
    { Histo & h = his_phi_topo[reg][topo];
      h.byVol = true; // PHYSICS
      h.init(0,maxDphi,nDphi, minKt,maxKt,nKt, minKt,maxKt,nKt,
             b+"rp_phi_"+regName[reg]+"_"+topos[topo]+".his"); }

    { Histo & h = his_phi[reg];
      h.byVol = true; // PHYSICS
      h.init(0,maxDphi,nDphi, minKt,maxKt,nKt, minKt,maxKt,nKt,
             b+"rp_phi_"+regName[reg]+".his"); }
  }
  fileReg.close();

  //
  for(int type = 0; type < nParts; type++)
  {
    int nmass = (type==pion ? nMass : (nMass*2)/5); // FIXME

    //
    for(int topo = 0; topo < nTopos; topo++)
    { Histo & h = his_mph_topo[type][topo];
      h.byVol = true; // PHYSICS
      h.init(0,maxMass,nmass, 0,maxDphi,nDphi, 
                              minKt,maxKt,nKt, minKt,maxKt,nKt,
             b+"mph_"+parts[type]+"_"+topos[topo]+".his");
             // "");
    }

    //
    for(int topo = 0; topo < nTopos; topo++)
    { Histo & h = his_sph_topo[type][topo];
      h.byVol = true; // PHYSICS 
      h.init(minSha,maxSha,nSha, 0,maxDphi,nDphi, 
                                 minKt,maxKt,nKt, minKt,maxKt,nKt,
             ""); }

    //
    { Histo & h = his_mph[type];
      h.byVol = true; // PHYSICS
      h.init(0,maxMass,nmass, 0,maxDphi,nDphi, 
                              minKt,maxKt,nKt, minKt,maxKt,nKt,
             b+"mph_"+parts[type]+".his"); }

    //
    { Histo & h = his_sph[type];
      h.byVol = true; // PHYSICS 
      h.init(minSha,maxSha,nSha, 0,maxDphi,nDphi,
                              minKt,maxKt,nKt, minKt,maxKt,nKt,
             ""); }

    //
    { Histo & h = his_mas[type];
      h.byVol = true; // PHYSICS
      h.init(0,maxMass,nmass, 
             minKt,maxKt,nKt, minKt,maxKt,nKt,
             b+"mas_"+parts[type]+".his"); }

    //
    { Histo & h = his_mhi[type];
      h.byVol = true; // PHYSICS
      h.init(0,maxMass,nmass/5,
             minKt,maxKt,nKt, minKt,maxKt,nKt,
             b+"mhi_"+parts[type]+".his"); }

    //
    { Histo & h = his_sha[type];
      h.byVol = true; // PHYSICS
      h.init(minSha,maxSha,nSha,
             minKt,maxKt,nKt, minKt,maxKt,nKt,
             b+"sha_"+parts[type]+".his"); }
  }

  // nonresonant continuum
  for(int type = 0; type < nParts; type++)
  {
    // topo
    for(int topo = 0; topo < nTopos; topo++)
    {

    { Histo & h = his_hat_topo[type][topo];
      h.byVol = true; // PHYSICS
      h.init(-10,0,(type==pion ?  50 : 25),
             -10,0,(type==pion ?  50 : 25), ""); }
    }

    // result
    { Histo & h = his_hat[type];
      h.byVol = true; // PHYSICS
      h.init(-10,0,(type==pion ?  50 : 25),
             -10,0,(type==pion ?  50 : 25), b+"hat_"+parts[type]+".his"); }
  }

}

/*****************************************************************************/
Distributions::~Distributions()
{
  cerr << Helper::col(5) << " writing distributions.." << Helper::col();

  //
  for(int ireg = 0; ireg < nRegs; ireg++)
  {
    vector<int> ix(3);

    for(ix[0] = 0; ix[0] < his_phi[ireg].axes[0].bins; ix[0]++) // dphi
    for(ix[1] = 0; ix[1] < his_phi[ireg].axes[1].bins; ix[1]++) // p1t
    for(ix[2] = 0; ix[2] < his_phi[ireg].axes[2].bins; ix[2]++) // p2t
    {
      average(his_phi_topo[ireg],ix, his_phi[ireg]); // over topos

      if(ix[1] < ix[2]) // p1t <-> p2t
      {
        for(int topo = 0; topo < nTopos; topo++)
        symmetrize(his_phi_topo[ireg][topo], ix);

        symmetrize(his_phi[ireg], ix);
      }
    }
  }

  ////////////////////////////////////////////////////////////////////
  for(int type = 0; type < nParts; type++)
  {
    vector<int> ix(4);
    for(ix[0] = 0; ix[0] < his_mph[type].axes[0].bins; ix[0]++) // m
    for(ix[1] = 0; ix[1] < his_mph[type].axes[1].bins; ix[1]++) // dphi
    for(ix[2] = 0; ix[2] < his_mph[type].axes[2].bins; ix[2]++) // p1t
    for(ix[3] = 0; ix[3] < his_mph[type].axes[3].bins; ix[3]++) // p2t
    {
      average(his_mph_topo[type],ix, his_mph[type]); // over topos

      if(ix[2] < ix[3]) // p1t <-> p2t
      {
        for(int topo = 0; topo < nTopos; topo++)
        symmetrize(his_mph_topo[type][topo], ix);

        symmetrize(his_mph[type], ix);
      }
    }

    // and fill his_mas
    for(ix[0] = 0; ix[0] < his_mph[type].axes[0].bins; ix[0]++) // m
    for(ix[2] = 0; ix[2] < his_mph[type].axes[2].bins; ix[2]++) // p1t
    for(ix[3] = 0; ix[3] < his_mph[type].axes[3].bins; ix[3]++) // p2t
    {
      float va=0, s2=0;
      for(ix[1] = 0; ix[1] < his_mph[type].axes[1].bins; ix[1]++) // dphi
      {
        va += his_mph[type].get(ix).first;
        s2 += his_mph[type].get(ix).second;
      }

      {
      vector<int> iy = {ix[0],ix[2],ix[3]};
      his_mas[type].get(iy) = pair<float,float>(va,s2);
      }

      { // hi
      vector<int> iy = {ix[0]/5,ix[2],ix[3]};
      his_mhi[type].get(iy).first  += va;
      his_mhi[type].get(iy).second += s2;
      }
    }

    //
    for(ix[0] = 0; ix[0] < his_sph[type].axes[0].bins; ix[0]++) // sha
    for(ix[1] = 0; ix[1] < his_sph[type].axes[1].bins; ix[1]++) // dphi
    for(ix[2] = 0; ix[2] < his_sph[type].axes[2].bins; ix[2]++) // p1t
    for(ix[3] = 0; ix[3] < his_sph[type].axes[3].bins; ix[3]++) // p2t
    {
      average(his_sph_topo[type],ix, his_sph[type]); // over topos

      if(ix[2] < ix[3]) // p1t <-> p2t
      {
        for(int topo = 0; topo < nTopos; topo++)
        symmetrize(his_sph_topo[type][topo], ix);

        symmetrize(his_sph[type], ix);
      }
    }

    // and fill his_sha
    for(ix[0] = 0; ix[0] < his_sph[type].axes[0].bins; ix[0]++) // m
    for(ix[2] = 0; ix[2] < his_sph[type].axes[2].bins; ix[2]++) // p1t
    for(ix[3] = 0; ix[3] < his_sph[type].axes[3].bins; ix[3]++) // p2t
    {
      float va=0, s2=0;
      for(ix[1] = 0; ix[1] < his_sph[type].axes[1].bins; ix[1]++) // dphi
      {
        va += his_sph[type].get(ix).first;
        s2 += his_sph[type].get(ix).second;
      }

      {
      vector<int> iy = {ix[0],ix[2],ix[3]};
      his_sha[type].get(iy) = pair<float,float>(va,s2);
      }
    }
  }

  // virtual
  for(int type = 0; type < nParts; type++)
  {
    { vector<int> ix(3);
      for(ix[0] = 0; ix[0] < his_hat[type].axes[0].bins; ix[0]++)
      for(ix[1] = 0; ix[1] < his_hat[type].axes[1].bins; ix[1]++)
        average(his_hat_topo[type],ix, his_hat[type]); }
  }

  cerr << " [done]" << endl;
}

/*****************************************************************************/
void Distributions::average(Histo his[], const vector<int> & ix,
                            Histo  & ave)
{
  vector<float> val,si2;
  for(int topo = 0; topo < nTopos; topo++)
  { 
    val.push_back(his[topo].get(ix).first );
    si2.push_back(his[topo].get(ix).second);
  }

  float va=0, s2=0; 
  for(int topo = 0; topo < nTopos; topo++)
  if(si2[topo] > 0)
  {
    va += val[topo]/si2[topo];
    s2 +=         1/si2[topo];
  }

  if(s2 >0) ave.get(ix) = pair<float,float>(va/s2,1/s2);
       else ave.get(ix) = pair<float,float>(0,0);
}

/*****************************************************************************/
// finally p1t < p2t
void Distributions::symmetrize(Histo & his, const vector<int> & ix)
{
  const int n = ix.size();
  vector<int> iy = ix; iy[n-2] = ix[n-1];
                       iy[n-1] = ix[n-2];

  vector<float> val,si2;

  val.push_back(his.get(ix).first );
  si2.push_back(his.get(ix).second);

  val.push_back(his.get(iy).first );
  si2.push_back(his.get(iy).second);

  float va=0, s2=0;
  for(int i = 0; i < 2; i++)
  if(si2[i] > 0)
  {
    va += val[i]/si2[i];
    s2 +=      1/si2[i];
  }

  if(s2 > 0) his.get(ix) = pair<float,float>(2 * va/s2,2 * 1/s2);
        else his.get(ix) = pair<float,float>(0,0);
}

/*****************************************************************************/
void Distributions::collectForPhysics(int topo, int type,
  float dphi_tilde,
  const vector<float> & dphip1tp2t, float m,
  const Vector4 & q1, const Vector4 & q2,
  const Vector4 & p3, const Vector4 & p4, float weight)
{
  const float & dphi = dphip1tp2t[0];
  const float &  p1t = dphip1tp2t[1];
  const float &  p2t = dphip1tp2t[2];

  //
  his_dphis[type].fillw({dphi_tilde,dphi},weight);

  // collect dphi in mass regions
  for(int ireg = 0; ireg < nRegs; ireg++)
  if(type == regType[ireg])
  if(lowSign[ireg] < m && m < higSign[ireg])
    his_phi_topo[ireg][topo].fillw({dphi,p1t,p2t}, weight);

  // collect (m,dphi)
  his_mph_topo[type][topo].fillw({m,dphi,p1t,p2t}, weight);

  // virtual
  Vector4 d;
  if(fabs((q1-p3).mass2()) <
     fabs((q1-p4).mass2())) d = q1 - p3;
                       else d = q1 - p4;

  float that_a = (q1-p3).mass2();
  float that_b = (q1-p4).mass2();

  if( (type == pion && m > 1.80 && m < 2.20) || // PARAMETERS
      (type == kaon && m > 2.20 && m < 2.60) || // PARAMETERS
      (type == prot && m > 0.00 && m < 9.99) )  // PARAMETERS
  {
    his_hat_topo[type][topo].fillw({    that_a,that_b }, weight);
    his_sph_topo[type][topo].fillw({max(that_a,that_b),dphi,p1t,p2t}, weight);
  }
}

