#ifndef _NLLS_H_
#define _NLLS_H_

#include <unsupported/Eigen/NonLinearOptimization>
#include "../EigenTypes.h"

template<int NX=Eigen::Dynamic, int NY=Eigen::Dynamic>
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


  NLLSObject() : inputs_(InputsAtCompileTime), values_(ValuesAtCompileTime) {}
  NLLSObject(unsigned int inputs, unsigned int values)
      : inputs_(inputs), values_(values)
  {}
  virtual ~NLLSObject() {}


  unsigned int inputs() const { return inputs_; }
  unsigned int values() const { return values_; }

  // you should define that in the subclass :
  virtual int operator() (const InputType& x, ValueType& v) const = 0;

 protected:
  const unsigned int inputs_, values_;

};

inline std::string lmInfo2String(int info)
{
  switch(info)
  {
    case -2: return "Not started";
    case -1: return "Running";
    case 0: return "Improper input parameters";
    case 1: return "Relative reduction too small";
    case 2: return "Relative error too small";
    case 3: return "Relative error and reduction too small";
    case 4: return "Cosinus too small";
    case 5: return "Too many function evaluations";
    case 6: return "Ftol too small";
    case 7: return "Xtol too small"; // <-- you want to have this
    case 8: return "Gtol too small";
    case 9: return "User asked";
    default: return "Unknown info";
  }
  return "Shouldn't happen";
}

#endif
