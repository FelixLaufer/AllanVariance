#ifndef _NLLS_FUNCTOR_H_
#define _NLLS_FUNCTOR_H_

#include "EigenTypes.h"
#include <unsupported/Eigen/NonLinearOptimization>
#include <unsupported/Eigen/NumericalDiff>

template<int NX = Eigen::Dynamic, int NY = Eigen::Dynamic>
class NLLSFunctor
{
 public:
  enum
  {
    InputsAtCompileTime = NX,
    ValuesAtCompileTime = NY
  };

  typedef ScalarType Scalar;
  typedef Eigen::Matrix<ScalarType, InputsAtCompileTime, 1> InputType;
  typedef Eigen::Matrix<ScalarType, ValuesAtCompileTime, 1> ValueType;
  typedef Eigen::Matrix<ScalarType, ValuesAtCompileTime, InputsAtCompileTime> JacobianType;

  NLLSFunctor()
    : inputs_(InputsAtCompileTime), values_(ValuesAtCompileTime)
  {}
 
  NLLSFunctor(unsigned int inputs, unsigned int values)
    : inputs_(inputs)
    , values_(values)
  {}
 
  virtual ~NLLSFunctor()
  {}

  unsigned int inputs() const
  {
    return inputs_;
  }
 
  unsigned int values() const
  {
    return values_;
  }
 
  virtual int operator() (const InputType& x, ValueType& fvec) const = 0;
  virtual int df(const InputType& x, JacobianType& fjac) const
  {};
 
 protected:
  const unsigned int inputs_, values_;
};

#endif
