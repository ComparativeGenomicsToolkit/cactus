// Copyright 2014 Martin C. Frith

// Fast calculation of: base ^ exponent
// where base is a "double" and exponent is an "int".

// Algorithm: "exponentiation by squaring".

// Assumes we will calculate many different exponents for the same base.

#ifndef CBRC_INT_EXPONENTIATOR
#define CBRC_INT_EXPONENTIATOR

namespace cbrc {

class IntExponentiator {
public:
  void setBase(double b) {
    base = b;
    invBase = 1.0 / b;
  }

  double operator()(int exponent) const {
    unsigned n = exponent;
    double x;
    if (exponent >= 0) {
      x = base;
    } else {
      x = invBase;
      n = -n;
    }
    double y = 1.0;
    while (n) {
      if (n % 2) y *= x;
      n /= 2;
      x *= x;
    }
    return y;
  }

private:
  double base;
  double invBase;
};

}

#endif
