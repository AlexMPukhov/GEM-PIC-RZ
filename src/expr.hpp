// ============================================================================
//  expr.hpp -- mathematical expressions in the input deck ("function parser")
//
//  An expression is compiled once into a small stack program and evaluated
//  for given values of its variables, e.g.
//      density_function = "exp(-(r/w)^2) * (1 + 0.1*cos(z/(50 um)))"
//
//  Syntax
//    numbers         2, 1.5e-3, and numbers with a unit: 10 um, 3e18 cm^-3, 5 fs
//                    (converted to code units like every value of the deck)
//    variables       given by the context (x, y, z, r, zeta, ...), code units
//    constants       pi, and the entries of the deck block `constants`
//    operators       + - * / ^ (also **), unary -, comparisons < <= > >= == !=
//                    (1 or 0), && || !, and  cond ? a : b
//    functions       sin cos tan asin acos atan atan2 sinh cosh tanh exp log
//                    log10 sqrt abs floor ceil min max pow sign step (1 for x >= 0)
//                    if(cond, a, b)  gauss(x, s) = exp(-x^2/(2 s^2))
//                    rect(x, w) = 1 for |x| <= w/2
//  Errors are thrown as std::invalid_argument with the position in the text.
// ============================================================================
#pragma once
#include <map>
#include <string>
#include <vector>

#include "units.hpp"

namespace gem {

class Expr {
 public:
  Expr() = default;
  // compile `text` with the variables `vars` (evaluated in this order), the
  // named constants `consts` (code units) and the unit conversion `units`
  Expr(const std::string& text, const std::vector<std::string>& vars,
       const std::map<std::string, double>& consts, const Units& units);

  double operator()(const double* v) const;
  double operator()(std::initializer_list<double> v) const { return (*this)(v.begin()); }
  bool empty() const { return code.empty(); }
  const std::string& text() const { return src; }

  enum Op : int {
    PUSH, VAR, NEG, NOT, ADD, SUB, MUL, DIV, POW, LT, LE, GT, GE, EQ, NE, AND, OR, SEL,
    F1, F2  // one/two-argument functions: index in `arg`
  };
  struct Ins {
    Op op;
    int arg;
    double val;
  };

 private:
  std::string src;
  std::vector<Ins> code;
  int maxStack = 0;
};

}  // namespace gem
