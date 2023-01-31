#ifndef _ALLAN_VARIANCE_H_
#define _ALLAN_VARIANCE_H_

#include "EigenTypes.h"
#include "NLLSFunctor.h"

#include <sstream>

namespace AllanVariance
{
  Vector generateOctaveSequence(const size_t len);
  Vector generateLogSpacedSequence(const size_t len, const size_t numMax);

  struct Allan
  {
    Vector variance;   // Allan variance, unit [unit^2]
    Vector deviation;  // Allan deviation, unit [unit]
    Vector tau;        // time series, unit [s]
    size_t numSamples; // Number of data samples used
    ScalarType dt;     // Sampling period
  };

  Allan compute(const Vector& data, const ScalarType& dt, const Vector& averagingFactors);
  Allan compute(const Vector& data, const ScalarType& dt, const size_t numAveragingFactors = 0);

  struct NoiseAnalysis
  {
    Allan allan;
    enum class Method { Slope = 0, Regression = 1 };
    Method method;       // Analysis method slope vs. reression approach
    ScalarType sigmaQ;   // Quantization error,  unit [unit*s]
    ScalarType sigmaRW;  // Random walk or noise desnity (white noise), unit [unit/sqrt(Hz)]
    ScalarType sigmaB;   // Bias instability (pink noise), unit [unit]
    ScalarType sigmaRRW; // Rate random walks (red or Brownian noise), unit [unit*sqrt(Hz)]
    ScalarType sigmaRR;  // Rate ramp, unit [unit/s] 
  };

  NoiseAnalysis analyzeSlopeMethod(const Allan& allan);

  struct RegressionFunctor : public NLLSFunctor<>
  {
    RegressionFunctor(const Matrix& X, const Vector& y)
      : NLLSFunctor<>(X.cols(), y.size())
      , X_(X)
      , y_(y)
    {}
  
    int operator()(const Vector& x, Vector& fvec) const
    {
      fvec = y_.array() - (X_ * x.cwiseAbs2()).array().log10();
      return 0;
    }

    Matrix X_;
    Vector y_;
  };

  struct RegressionFunctorNumDiff : public Eigen::NumericalDiff<RegressionFunctor> {};

  NoiseAnalysis analyzeRegressionMethod(const Allan& allan);

  std::stringstream write(const NoiseAnalysis& result);
}

#endif
