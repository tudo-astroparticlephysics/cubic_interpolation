#pragma once

#include <cmath>
#include <ostream>

class Axis {
protected:
  float low;
  float stepsize;

  virtual void print(std::ostream &) const = 0;

public:
  Axis(float _low, float _stepsize) : low(_low), stepsize(_stepsize) {}

  virtual size_t node(float &x) {
    size_t n = std::floor((x - low) / stepsize);
    x -= low + n * stepsize;
    return n;
  };
  virtual float back_node(size_t n) { return low + n * stepsize; };
  virtual void transform(float &) noexcept = 0;
  virtual void back_transform(float &) noexcept = 0;

  friend std::ostream &operator<<(std::ostream &out, const Axis &);
};

struct LogAxis : public Axis {
  LogAxis(float _low, float _stepsize) : Axis(_low, _stepsize) {}
  inline void transform(float &x) noexcept final { x = std::log(x); }
  inline void back_transform(float &x) noexcept final { x = std::exp(x); }

protected:
  void print(std::ostream &) const;
};

struct LinAxis : public Axis {
  LinAxis(float _low, float _stepsize) : Axis(_low, _stepsize) {}
  inline void transform(float &) noexcept final { return; }
  inline void back_transform(float &) noexcept final { return; }

protected:
  void print(std::ostream &) const;
};

std::ostream &operator<<(std::ostream &out, const Axis &axis) {
  axis.print(out);
  return out;
}

void LogAxis::print(std::ostream &os) const {
  os << "LogAxis(low: " << low << ", stepsize: " << stepsize << ")";
}

void LinAxis::print(std::ostream &os) const {
  os << "LinAxis(low: " << low << ", stepsize: " << stepsize << ")";
}
