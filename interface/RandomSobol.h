#ifndef _RandomSobol_h_
#define _RandomSobol_h_

#include <vector>
#include <map>

class RandomSobol
{
 public:
  RandomSobol(const unsigned int N, const unsigned int D);
  ~RandomSobol() {};

  void reset();
  std::vector<double> get();

 private:
  const unsigned int N;
  const unsigned int D;

  std::vector<unsigned int> C;
  std::vector<std::map<unsigned int, unsigned int>> V; // [D][L+1];

  std::vector<unsigned int> X;

  int i;
};

#endif
