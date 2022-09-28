#ifndef _NLLS_H_
#define _NLLS_H_

#include "EigenTypes.h"
#include <unsupported/Eigen/NonLinearOptimization>

template<int NX = Eigen::Dynamic, int NY = Eigen::Dynamic>
class NLLSObject
{
 public:
  enum
  {
    InputsAtCompileTime = NX,
    ValuesAtCompileTime = NY
  };

  typedef ScalarType Scalar;
  typedef Eigen::Matrix<ScalarType,InputsAtCompileTime,1> InputType;
  typedef Eigen::Matrix<ScalarType,ValuesAtCompileTime,1> ValueType;
  typedef Eigen::Matrix<ScalarType,ValuesAtCompileTime,InputsAtCompileTime> JacobianType;

  NLLSObject()
    : inputs_(InputsAtCompileTime), values_(ValuesAtCompileTime)
  {}
 
  NLLSObject(unsigned int inputs, unsigned int values)
    : inputs_(inputs), values_(values)
  {}
 
  virtual ~NLLSObject() {}

  unsigned int inputs() const { return inputs_; }
  unsigned int values() const { return values_; }
 
  virtual int operator() (const InputType& x, ValueType& v) const = 0;

 protected:
  const unsigned int inputs_, values_;
};

#endif
