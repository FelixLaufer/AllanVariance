#ifndef _ALLAN_VARIANCE_H_
#define _ALLAN_VARIANCE_H_

#include "math/EigenTypes.h"
#include "math/optim/EigenNLLS.h"

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

  struct RegressionFunctor : public NLLSObject<>
  {
    int operator()(const Vector& x, Vector& fvec) const
    {
      fvec(0) = (y.array() - (X * x.cwiseAbs2()).array().log10()).cwiseAbs2().sum();
      return 0;
    }

    unsigned int inputs() const { return 5; };
    unsigned int values() const { return 1; };

    Matrix X;
    Vector y;
  };

  struct RegressionFunctorNumDiff : public Eigen::NumericalDiff<RegressionFunctor> {};

  NoiseAnalysis analyzeRegressionMethod(const Allan& allan);

  std::stringstream write(const NoiseAnalysis& result);
}

#endif