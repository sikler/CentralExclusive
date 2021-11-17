// https://web.maths.unsw.edu.au/~fkuo/sobol/
#include "../interface/RandomSobol.h"

#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

/*****************************************************************************/
RandomSobol::RandomSobol(const unsigned int N, const unsigned int D)
                        : N(N), D(D)
{
//  cerr << " initializing RandomSobol's sequence (" << N << "," << D << ")..";

  // L = max number of bits needed 
  unsigned int L = (unsigned int)ceil(log((double)N)/log(2.0));

  // C[i] = index from the right of the first zero bit of i
  C.resize(N);

  C[0] = 1;
  for (unsigned int i=1;i<N;i++) {
    C[i] = 1;
    unsigned int value = i;
    while (value & 1) {
      value >>= 1;
      C[i]++;
    }
  }

  // ----- Compute the first dimension -----
  V.resize(D);

  // Compute direction numbers V[1] to V[L], scaled by pow(2,32)
  for (unsigned i=1;i<=L;i++) V[0][i] = 1 << (32-i); // all m's = 1

  // ----- Compute the remaining dimensions -----
  ifstream file("../pars/new-joe-kuo-6.21201");

  char buffer[1000];
  file.getline(buffer,1000,'\n');

  for(unsigned j=1;j<D;j++)
  {
    // Read in parameters from file 
    unsigned int d, s, a;
    file >> d >> s >> a;

    vector<unsigned int> m(s+1);
    for (unsigned int i=1;i<=s;i++) file >> m[i];

    // Compute direction numbers V[1] to V[L], scaled by pow(2,32)
    if (L <= s) {
      for (unsigned int i=1;i<=L;i++) V[j][i] = m[i] << (32-i); 
    }
    else {
      for (unsigned int i=1;i<=s;i++) V[j][i] = m[i] << (32-i); 
      for (unsigned int i=s+1;i<=L;i++) {
        V[j][i] = V[j][i-s] ^ (V[j][i-s] >> s); 
        for (unsigned int k=1;k<=s-1;k++)
          V[j][i] ^= (((a >> (s-1-k)) & 1) * V[j][i-k]); 
      }
    }
  }

  file.close();

  i = 0;

  X.resize(D,0);

//  cerr << " [done]" << endl;
}

/*****************************************************************************/
void RandomSobol::reset()
{
  i = 0;
  fill(X.begin(), X.end(), 0);
}

/*****************************************************************************/
vector<double> RandomSobol::get()
{
  vector<double> POINTS(D);

  for(unsigned int j=0;j<D;j++)
  {
    X[j] = X[j] ^ V[j][C[i-1]];
    POINTS[j] = (double)X[j]/pow(2.0,32); // *** the actual points
  }

  i++;

  return POINTS;
}

