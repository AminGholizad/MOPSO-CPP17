#ifndef PSO_H
#define PSO_H
#include <cmath>
#include "particle.hpp"
#include "repository.hpp"
namespace mopso {
  template<ull N,ull O,ull S=100,ull R=100>
  std::pair<Repository<N,O,S>,std::array<Particle<N,O,S>,S>> mopso(const vars<N>& lower_bound,const vars<N>& upper_bound,const Problem<N,O>& problem,const size_t max_iter=1000,const std::array<double,2>& c={0.2,0.2},const std::array<double,2>& iw={0.1,0.01},const double mu=0.1){
    auto w = [&](size_t it){ return ((max_iter - it) - (iw[0] - iw[1]))/max_iter + iw[1];};
    auto pm = [&](size_t it){ return std::pow(1-it/(max_iter-1.),1/mu);};
    std::array<Particle<N,O,S>,S> swarm;
    for (size_t i = 0; i < S; i++)
      swarm[i]=Particle<N,O,S>(lower_bound,upper_bound,problem);
    Repository<N,O,S> rep(swarm,R,7,0.1,2,2);
    for (size_t i = 0; i < max_iter; i++) {
      Particle<N,O,S> gBest = rep.SelectLeader();
      auto wc = w(i);
      auto pc = pm(i);
      for (size_t j = 0; j < S; j++) {
        swarm[j].update(gBest,problem,wc,c,pc);
      }
      rep.update(swarm);
    }
    return std::make_pair(rep,swarm);
  }
} /* pso */
#endif //PSO_H
