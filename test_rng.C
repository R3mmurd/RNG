/* Algorithms to statistic tests for pseudo random number generators.

   This program contains:

  - Kolmogorov-Smirnov algorithm to check uniformity.
  - Breusch-Godfrey algorithm to test serial correlation.
  - Serial test for 2-dimension uniformity.
  - Streak test.
  - Gap test.
  - Poker test.

   Autor: Alejandro J. Mujica (aledrums en gmail.com)
*/

# include <cstdlib>
# include <cmath>
# include <cassert>

# include <iostream>
# include <sstream>
# include <iomanip>
# include <vector>
# include <algorithm>

# include <rng.H>

// Chi square table with alpha = 0.05
double chi_squared(const unsigned int & v)
{
  switch(v)
    {
    case 1:          return 3.84;
    case 2:          return 5.99;
    case 3:          return 7.81;
    case 4:          return 9.49;
    case 5:          return 11.07;
    case 6:          return 12.59;
    case 7:          return 14.07;
    case 8:          return 15.51;
    case 9:          return 16.92;
    case 10:         return 18.31;
    case 11:         return 19.68;
    case 12:         return 21.03;
    case 13:         return 22.36;
    case 14:         return 23.68;
    case 15:         return 25.0;
    case 16:         return 26.3;
    case 17:         return 27.59;
    case 18:         return 28.87;
    case 19:         return 30.14;
    case 20:         return 31.41;
    case 21:         return 32.67;
    case 22:         return 33.92;
    case 23:         return 35.17;
    case 24:         return 36.42;
    case 25:         return 37.65;
    case 26:         return 38.89;
    case 27:         return 40.11;
    case 28:         return 41.34;
    case 29:         return 42.56;
    case 30:         return 43.77;
    case 31 ... 40:  return 55.76;
    case 41 ... 50:  return 67.5;
    case 51 ... 60:  return 79.08;
    case 61 ... 70:  return 90.53;
    case 71 ... 80:  return 101.88;
    case 81 ... 90:  return 113.15;
    case 91 ... 100: return 124.34;
    default:         return 0.0;
    }
}

// Kolmogorov-Smirnov test table
double dtab(const size_t & n)
{
  switch(n)
    {
    case 1 ... 5:   return 0.565;
    case 6:         return 0.521;
    case 7:         return 0.486;
    case 8:         return 0.457;
    case 9:         return 0.432;
    case 10:        return 0.410;
    case 11:        return 0.391;
    case 12:        return 0.375;
    case 13:        return 0.361;
    case 14:        return 0.349;
    case 15:        return 0.338;
    case 16:        return 0.328;
    case 17:        return 0.318;
    case 18:        return 0.309;
    case 19:        return 0.301;
    case 20:        return 0.294;
    case 21 ... 25: return 0.27;
    case 26 ... 30: return 0.24;
    case 31 ... 35: return 0.23;
    default:        return 1.36 / std::sqrt(n);
    }
}

// Uniformity test
void Kolmogorov_Smirnov_test(std::vector<double> & a)
{
  std::sort(a.begin(), a.end());

  const size_t & n = a.size();

  double inv_n = 1.0 / double(n);

  double dplus = std::abs(inv_n - a[0]);

  double dminus = a[0]; // abs(a[0] - 0 / n)

  for (size_t i = 2; i < n; ++i)
    {
      double aux = std::abs(i * inv_n - a[i - 1]);

      if (aux > dplus)
        dplus = aux;

      aux = std::abs(a[i - 1] - (i - 1) * inv_n);

      if (aux > dminus)
        dminus = aux;
    }

  double dt = dtab(n);

  const double & dc = std::max(dplus, dminus);

  std::cout << "----- UNIFORMITY TEST (KS) -----\n\n"
            << "Calculated D  Tabulated D (0.05)\n"
            << "------------  ------------------\n"
            << std::setprecision(4) << std::setw(14) << dc
            << std::setw(20) << dt;

  if (dc < dt)
    std::cout << "PASSED\n";
  else
    std::cout << "FAILED\n";
}

// Serial correlation test
void Breusch_Godfrey_test(std::vector<double> & a)
{
  const size_t & n = a.size();

  const size_t max_k = 10;

  std::cout << "-------------- SERIAL CORRELATION TEST ---------------\n\n"
            << "Distance   Autocovariance     Confidence Interval 95%\n"
            << "    k           Rk          Lower bound   Upper bound\n"
            << "---------  --------------   -----------   -----------\n";

  for (size_t k = 1; k <= max_k; ++k)
    {
      double sum = 0;

      for (size_t i = 0; i < n - k; ++i)
        sum += (a[i] - 0.5) * (a[i + k] - 0.5);

      double rk = sum / (n - k);

      // Intervalo de confianza
      double v = 1.958 / (12 * std::sqrt(n - k));
      double li = rk - v;
      double ls = rk + v;

      std::cout << std::fixed << std::setprecision(7)
                << std::setw(11) << k << std::setw(17) << rk
                << std::setw(14) << li << std::setw(14) << ls;

      if (li < 0 and ls > 0)
        std::cout << " PASSED\n";
      else
        std::cout << " FAILED\n";
    }
}

// 2-dimensional uniformity test.
void serial_test_2D(std::vector<double> & a)
{
  constexpr size_t cell_size = 5;

  unsigned int matrix[cell_size][cell_size];

  for (size_t i = 0; i < cell_size; ++i)
    for (size_t j = 0; j < cell_size; ++j)
      matrix[i][j] = 0;

  size_t num_cells = cell_size * cell_size;

  for (size_t i = 0; i < a.size(); i += 2)
    {
      double x = a[i];
      double y = a[i + 1];

      size_t m_i = y * cell_size;
      size_t m_j = x * cell_size;

      ++matrix[m_i][m_j];
    }

  double expected = double(a.size()) / (2.0 * double(num_cells));

  double kc = 0.0;

  for (size_t i = 0; i < cell_size; ++i)
    for (size_t j = 0; j < cell_size; ++j)
      {
        double diff = double(matrix[i][j]) - expected;
        kc += ((diff * diff) / expected);
      }

  double kt = chi_squared(num_cells - 1);

  std::cout << "-- 2D UNIFORMITY TEST --\n\n"
            << "Calculated K  Tabulated K (0.05)\n"
            << "------------  ------------------\n"
            << std::setprecision(4) << std::setw(14) << kc
            << std::setw(19) << kt;


  if (kc < kt)
    std::cout << "PASSED\n";
  else
    std::cout << "FAILED\n";
}

void streak_test(std::vector<double> & a)
{
  const double Z_005 = 1.645;

  std::vector<char> streaks;

  const size_t & n = a.size();

  for (size_t i = 0; i < n - 1; ++i)
    {
 
      if (a[i] < a[i + 1])
        streaks.push_back('-');
      else
        streaks.push_back('+');
    }

  size_t num_streaks = 1;

  char curr_sign = streaks[0];

  for (const char & sign : streaks)
    {
      if (sign == curr_sign)
        continue;
        
      ++num_streaks;
      curr_sign = sign;
    }

  double mu = (2.0 * double(n) - 1) / 3.0;

  double sigma = (16.0 * double(n) - 29) / 90;

  double Z = (num_streaks - mu) / sigma;

  std::cout << "------------------ STREAK TEST -------------------\n\n"
            << "   a        mu       sigma        Z        Z(0.05)\n"
            << "-------   -------  ---------  ---------  ---------\n"
            << std::setprecision(4) << std::left
            << std::setw(10) << num_streaks << std::setw(9) << mu
            << std::setw(11) << sigma << std::setw(11) << Z << std::setw(9)
            << Z_005;

  if (Z > -Z_005 and Z < Z_005)
    std::cout << "PASSED\n";
  else
    std::cout << "FAILED\n";
}


void gap_test(std::vector<double> & a)
{
  const size_t num_digits = 10;

  std::vector<int> aux;

  for (const double & n : a)
    aux.push_back(int(n * num_digits));

  const size_t gap_interval = 4;

  size_t max_gap_size = aux.size() - 2;

  size_t num_gaps = aux.size() - num_digits;

  // Tamaño máximo del arreglo de frecuencias
  size_t freq_size = (max_gap_size / gap_interval) +
                     (max_gap_size % gap_interval == 0 ? 0 : 1);

  std::vector<int> gap_frequencies(freq_size);

  size_t max_pos = 0;

  for (size_t i = 0; i < num_digits; ++i)
    {
      // Ubicar el primer valor igual a i.

      size_t j;

      for (j = 0; j < aux.size() and aux[j] != i; ++j);

      // Tomar tamaños de brechas
      size_t gap_length = 0;

      ++j;

      for ( ; j < aux.size(); ++j)
        {
          if (aux[j] == i)
            {
              size_t freq_pos = gap_length / gap_interval;
              max_pos = std::max(freq_pos, max_pos);
              ++gap_frequencies[freq_pos];
              gap_length = 0;
            }
          else
            ++gap_length;
        }
    }

  const double D_005 = 0.136;

  std::cout << "------------------------------ GAP TEST -----------------------"
	    << "--------\n"
            << "                                     Cumulative\n"
	    << "                                     Relative\n"
            << "                         Relative    Frequency\n"
            << "Gap length  Frequency    Frequency      S(x)      F(x)    "
	    << "|F(x) - S(x)|\n"
            << "----------  ----------   ----------  ----------  ------   "
	    << "-------------\n";

  double last_rel_freq = 0.0;

  double max_diff = 0.0;

  for (size_t i = 0; i <= max_pos; ++i)
    {
      double relative_freq = double(gap_frequencies[i]) / double(num_gaps);

      size_t lr = i * gap_interval;

      size_t ll = i * gap_interval + gap_interval - 1;

      double sx = relative_freq + last_rel_freq;

      last_rel_freq = sx;

      double fx = 1 - std::pow(0.9, ll + 1);

      double fxsx = std::abs(fx - sx);

      std::stringstream sstr;

      sstr << lr << '-' << ll;

      std::cout << std::left << std::setw(12) << sstr.str()
                << std::setw(13) << gap_frequencies[i]
                << std::setw(12) << relative_freq << std::setw(12)
                << sx << std::setw(9) << fx << std::setw(13) << fxsx;

      std::cout << '\n';

      max_diff = std::max(fxsx, max_diff);
    }

    std::cout << "Maximum difference: " << max_diff << '\n';
    std::cout << "Critial value: " << D_005 << '\n';
    std::cout << "Result: ";

    if (max_diff < D_005)
      std::cout << "PASSED\n";
    else
      std::cout << "FAILED\n";
}

void poker_test(std::vector<double> & a)
{
  const size_t & n = a.size();

  size_t none = 0;
  size_t all  = 0;
  size_t pair = 0;

  const double expected_none = n * 0.72;
  const double expected_all  = n * 0.01;
  const double expected_pair = n * 0.27;

  for (const double & sample : a)
    {
      std::stringstream sstr;
      sstr << std::fixed << std::setprecision(3) << sample;

      std::string str = sstr.str();

      assert(str.size() == 5);

      if (str[3] == str[2] and str[4] == str[3])
        ++all;
      else if (str[3] == str[2] or str[4] == str[3] or str[4] == str[2])
        ++pair;
      else
        ++none;
    }

    double op_none = (none - expected_none) *
                     (none - expected_none) / expected_none;

    double op_all  = (all - expected_all) *
                     (all - expected_all) / expected_all;

    double op_pair = (pair - expected_pair) *
                     (pair - expected_pair) / expected_pair;

    size_t O = none + all + pair;

    double E = expected_none + expected_all + expected_pair;

    double T = op_none + op_all + op_pair;

    std::cout << "----------------------- POKER TEST -----------------------\n"
              << "              Observed         Expected        (Oi - Ei)^2\n"
              << "             Frequencies      Frequencies      -----------\n"
              << "Case              Oi              Ei               Ei     \n"
              << "----------------------------------------------------------\n"
              << "Distincts    " << std::setw(17) << none
                                 << std::setw(17) << expected_none
                                 << op_none << '\n'
              << "Equals       " << std::setw(17) << all
                                 << std::setw(17) << expected_all
                                 << op_all << '\n'
              << "A pair       " << std::setw(17) << pair
                                 << std::setw(17) << expected_pair
                                 << op_pair << '\n'
              << "----------------------------------------------------------\n"
              << "             " << std::setw(17) << O
                                 << std::setw(17) << E
                                 << T << '\n';

  const double X = 5.99;

  std::cout << "Chi-square(0.05,2) = " << X << '\n';

  if (T < X)
    std::cout << "PASSED\n";
  else
    std::cout << "FAILED\n";
}

int main(int argc, char * argv[])
{
  // Cantidad de pruebas que se le harán.
  unsigned long num_samples = argc < 2 ? 100 : std::stoul(argv[1]);

  long long seed = argc < 3 ? 1 : std::stoll(argv[1]);

  using Rng = LCG<long long>;

  seed %= Rng::max();

  // Generador sometido a pruebas. (semilla por omisión = 1)
  Rng lcg(seed);

  std::cout << "\n\nSample size: " << num_samples << "\n\n";

  std::vector<double> samples;

  for (unsigned long j = 0; j < num_samples; ++j)
    samples.push_back(lcg());

  streak_test(samples);

  std::cout << '\n';

  gap_test(samples);

  std::cout << '\n';

  Breusch_Godfrey_test(samples);

  std::cout << '\n';

  serial_test_2D(samples);

  std::cout << '\n';

  Kolmogorov_Smirnov_test(samples);

  std::cout << '\n';

  poker_test(samples);

  return 0;
}

