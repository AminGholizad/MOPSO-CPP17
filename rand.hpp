#ifndef RAND
#define RAND
#include <random>
#include <algorithm>
#include <chrono>
namespace rnd{
  inline std::mt19937& Generator(){
    static std::random_device rd;
    static unsigned seed = (rd.entropy()==0)?std::chrono::system_clock::now().time_since_epoch().count():rd();
    static std::mt19937 gen(seed);
    return gen;
  }
  template <typename T,typename RandomGenerator>
  inline T unifrnd(T a,T b,RandomGenerator& g){
    if constexpr(std::is_integral_v<T>){
      std::uniform_int_distribution<T> dis(a, b);
      return dis(g);
    }
    else {
      std::uniform_real_distribution<T> dis(a, b);
      return dis(g);
    }
  }
  template <typename T>
  inline T unifrnd(T a,T b){
    return unifrnd<T>(a,b,Generator());
  }
  inline double rand(){
    return unifrnd<double>(0.,1.);
  }
  template<typename Iter, typename RandomGenerator>
  inline Iter select_randomly(Iter start, Iter end, RandomGenerator& g){
      std::uniform_int_distribution<int> dis(0, std::distance(start, end) - 1);
      std::advance(start, dis(g));
      return start;
  }
  template<typename Iter>
  inline Iter select_randomly(Iter start, Iter end){
    return select_randomly(start,end,Generator());
  }
}
#endif
