# ifndef RNG_H
# define RNG_H

# include <limits>

/** Template class that implements a multiplicative linear congruential
 *  generator for pseudo-random numbers.
 *
 *  @tparam T The data type for seed. 
 *  @tparam A The multiplier
 */
template <typename T, T A = 16807>
class LCG
{
  static_assert(std::is_integral<T>::value,
                "template argument not an integral type");

  static constexpr T M = std::numeric_limits<T>::max();
  static constexpr T Q = M / A;
  static constexpr T R = M % A;

  T seed;

public:
  LCG(const T & _seed = T(time(0) % M))
    : seed(_seed)
  {
    // Empty
  }

  /// Return the max value generated.
  static constexpr T max()
  {
    return LCG::M;
  }

  /// The a seed.
  void set(const T & _seed)
  {
    seed = T(_seed);
  }

  /// Return the next generated seed.
  const T & get()
  {
    seed = A * (seed % Q) - R * (seed / Q);

    if (seed < 0)
      seed += M;

    return seed;
  }

  /// Return a pseudo-random number in [0, 1).
  double operator () ()
  {
    return double(get()) / double(M);
  }
};

# endif // RNG_H

