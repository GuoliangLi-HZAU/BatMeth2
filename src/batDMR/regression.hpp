
// Declares the Regression class which implements beta-binomial regression.

#ifndef REGRESSION_HPP_
#define REGRESSION_HPP_

#include <vector>
#include <string>

struct Design {
  std::vector<std::string> factor_names;
  std::vector<std::string> sample_names;
  std::vector<std::vector<double> > matrix;
};

std::istream& operator>> (std::istream &is, Design &design);
std::ostream& operator<< (std::ostream &os, const Design &design);
void remove_factor(Design &design, size_t factor);

struct SiteProportions {
  std::string chrom;
  size_t position;
  std::string strand;
  std::string context;
  std::vector<size_t> meth;
  std::vector<size_t> total;
  std::string name;
};

std::istream&
operator>>(std::istream &table_encoding, SiteProportions &props);

struct Regression {
  Design design;
  SiteProportions props;
  double max_loglik;
};

bool fit(Regression &r,
          std::vector<double> initial_parameters = std::vector<double>());

#endif // REGRESSION_HPP_
