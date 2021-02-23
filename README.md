# cubic_interpolation [![Build Status](https://api.travis-ci.com/MaxSac/cubic_interpolation.svg?branch=main)](https://travis-ci.com/github/MaxSac/cubic_interpolation)

A leightweight interpolation library based on boost and eigen.
It's provide the utilities to handle tables which are runtime intensive to build
and reduces interpolation failures to axis transformation.

If runtime intensive interpolation tables are created over several orders of magnitude,
this is an ideal application example.

```
#include <CubicInterpolation/Axis.h>
#include <CubicInterpolation/CubicSplines.h>
#include <CubicInterpolation/Interpolant.h>

using namespace cubic_splines;

auto lower_lim = 1.f;
auto upper_lim = 1e14.f;
auto nodes = 100u;

auto def = CubicSplines<double>::Definition();     // container of interpolation settings
def.f = func;                                      // runtime intensive function to approximate
def.axis = std::make_unique<LinAxis<double>>(
    lower_lim, upper_lim, nodes);                  // axis over many order of magnitudes

// interpolant which stores the tables at relocateable place to save time by
// second call.
auto inter = Interpolant<CubicSplines<double>>(std::move(def), TABLS_PATH, TABLES_NAM);

// evaluation of interpolation doing all neccessary retransformations
auto res = inter.evaluate(point);
```

More information can be found in the documentation which could be build with the
flag `BUILD_DOCUMENTATION` with a sphinx and a doxygen target. 


