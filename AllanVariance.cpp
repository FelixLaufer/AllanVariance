#include "AllanVariance.h"

#include "math/Constants.h"

#include <numeric>
#include <vector>

Vector AllanVariance::generateOctaveSequence(const size_t len)
{
  std::vector<ScalarType> ret;
  const size_t max = std::pow(2, std::floor(std::log2(static_cast<ScalarType>(len - 1) / 2)));
  for (size_t n = 0, m = 0; m <= max; ++n)
    ret.push_back(m = std::pow(2, n));
  return Vector(Eigen::Map<Vector>(ret.data(), ret.size()));
}

Vector AllanVariance::generateLogSpacedSequence(const size_t len, const size_t numMax)
{
  const size_t max = std::pow(2, std::floor(std::log2(static_cast<ScalarType>(len - 1) / 2)));
  const ScalarType step = std::log10(max) / (numMax - 1);
  Vector ret = Vector::Ones(numMax) * step;
  ret(0) = 0;
  std::partial_sum(ret.data(), ret.data() + ret.size(), ret.data(), std::plus<ScalarType>());
  ret.array() = (10 * Vector::Ones(numMax)).array().pow(ret.array()).ceil();
  return Vector(Eigen::Map<Vector>(ret.data(), std::unique(ret.data(), ret.data() + numMax) - ret.data()));
}

AllanVariance::Allan AllanVariance::compute(const Vector& data, const ScalarType& dt, const Vector& avrFcts)
{
  const size_t len = data.size();
  const Vector tau = dt * avrFcts;
  Vector theta = Vector::Zero(data.size());
  std::partial_sum(data.data(), data.data() + data.size(), theta.data(), std::plus<ScalarType>());
  theta *= dt;

  Vector avar = Vector::Zero(avrFcts.size());
  #pragma omp parallel for
  for (int i = 0; i < avrFcts.size(); ++i)
    avar(i) = (theta.segment(2 * avrFcts(i), len - 2 * avrFcts(i)) - 2 * theta.segment(avrFcts(i), len - 2 * avrFcts(i)) + theta.segment(0, len - 2 * avrFcts(i))).array().pow(2).sum();
  avar.array() = avar.array().cwiseQuotient(2 * tau.array().pow(2) * ((Vector::Ones(avrFcts.size()).array() * len) - 2 * avrFcts.array()));
  
  return { avar, avar.cwiseSqrt(), tau, len, dt };
}

AllanVariance::Allan AllanVariance::compute(const Vector& data, const ScalarType& dt, const size_t numAveragingFactors)
{
  return compute(data, dt, (numAveragingFactors == 0 ? generateOctaveSequence(data.size()) : generateLogSpacedSequence(data.size(), numAveragingFactors)));
}

AllanVariance::NoiseAnalysis AllanVariance::analyzeSlopeMethod(const Allan& allan)
{
  const Vector logADev = allan.deviation.array().log10();
  const Vector logTau = allan.tau.array().log10();
  const Vector dLogADev = logADev.segment(1, logADev.size() - 1) - logADev.segment(0, logADev.size() - 1);
  const Vector dLogTau = logTau.segment(1, logTau.size() - 1) - logTau.segment(0, logTau.size() - 1);

  const auto intercept = [&logADev, &dLogADev, &logTau, &dLogTau](const ScalarType& slope)
  {
    Eigen::Index minIdx;
    ((dLogADev.cwiseQuotient(dLogTau).array() - slope).cwiseAbs()).minCoeff(&minIdx);
    return logADev(minIdx) - slope * logTau(minIdx);
  };

  const auto sigmaFromSlopeAndTau = [&intercept](const ScalarType& slope, const ScalarType& tau)
  {
    return std::pow(10, slope * std::log10(tau) + intercept(slope));
  };
  
  const ScalarType sigmaQ = sigmaFromSlopeAndTau(-1, std::sqrt(3));
  const ScalarType sigmaRW = sigmaFromSlopeAndTau(-0.5, 1);
  const ScalarType sigmaB = std::pow(10, intercept(0) - std::log10(std::sqrt(2 * std::log(2) / Pi)));
  const ScalarType sigmaRRW = sigmaFromSlopeAndTau(0.5, 3);
  const ScalarType sigmaRR = sigmaFromSlopeAndTau(1, std::sqrt(2));

  return { allan, NoiseAnalysis::Method::Slope, sigmaQ, sigmaRW, sigmaB, sigmaRRW , sigmaRR };
}

AllanVariance::NoiseAnalysis AllanVariance::analyzeRegressionMethod(const Allan& allan)
{
  const size_t len = std::min(allan.deviation.size(), allan.tau.size());
  const Vector& aDev = allan.deviation;
  const Vector& aVar = allan.variance;
  const Vector& tau = allan.tau;

  // Build regressor matrix
  Matrix X(len, 5);
  X << 3 / tau.cwiseAbs2().array(), 1 / tau.array(), Vector::Ones(len), tau.array() / 3, tau.cwiseAbs2().array() / 2;
  
  // Get initual estimate from ridge regression
  const Matrix& X1 = X.cwiseSqrt();
  const ScalarType lambda = 5e-3;
  const Eigen::JacobiSVD<Matrix>& svd = X1.jacobiSvd(Eigen::ComputeFullU | Eigen::ComputeFullV);
  const Vector& s = svd.singularValues();
  const Eigen::Index r = s.rows();
  const Vector& x0 = svd.matrixV().leftCols(r) * s.cwiseQuotient((s.array().square() + lambda).matrix()).asDiagonal() * svd.matrixU().transpose().topRows(r) * aDev;

  // Get final estimate from non-linear least-squares optimization
  RegressionFunctorNumDiff func;
  const Matrix& X2 = X;
  func.X = X2;
  func.y = aVar.array().log10();
  Eigen::LevenbergMarquardt<RegressionFunctorNumDiff, ScalarType> lm(func);
  Vector x1 = x0;
  const int ret = lm.minimize(x1);

  // Extract sigmas
  Vector sigmas = x1.cwiseAbs();
  sigmas(2) = (X2 * sigmas.cwiseAbs2()).cwiseSqrt().minCoeff() / std::sqrt(2 * std::log2(2) / Pi);

  return { allan, NoiseAnalysis::Method::Regression, sigmas(0), sigmas(1), sigmas(2), sigmas(3), sigmas(4) };
}

std::stringstream AllanVariance::write(const NoiseAnalysis& noiseAnalysis)
{
  std::stringstream ret;
  ret << "Allan variance noise analysis ";
  ret << "using " << (noiseAnalysis.method == NoiseAnalysis::Method::Slope ? "slope" : "regression") << " method ";
  ret << "with " << noiseAnalysis.allan.numSamples << " data samples and sampling period " << noiseAnalysis.allan.dt << " [s]\n";
  ret << "sigmaQ\t\tquantization error\t" << noiseAnalysis.sigmaQ << "\t[unit*s]\n";
  ret << "sigmaRW\t\trandom walk\t\t" << noiseAnalysis.sigmaRW << "\t[unit/sqrt(Hz)]\n";
  ret << "sigmaB\t\tbias instability\t" << noiseAnalysis.sigmaB << "\t[unit]\n";
  ret << "sigmaRRW\trate random walk\t" << noiseAnalysis.sigmaRRW << "\t[unit*sqrt(Hz)]\n";
  ret << "sigmaRR\t\tramp rate\t\t" << noiseAnalysis.sigmaRR << "\t[unit/s]\n";
  return ret;
}