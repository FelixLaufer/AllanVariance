# AllanVariance
Allan variance and noise analysis for computing quantization error, random walk, bias instability rate random walk and rate ramp coefficients based on
- a) the common slope method and
- b) a regression approach as introduced in the paper "A regression-based methodology to improve estimation of inertial sensor errors using Allan variance data" by Jurado et al. (https://www.researchgate.net/publication/330514910_A_regression-based_methodology_to_improve_estimation_of_inertial_sensor_errors_using_Allan_variance_data)

## Example

const Vector data = Vector::Random(10000);
const auto& allanRes = AllanVariance::compute(data, 0.01);
const auto& slopeRes = AllanVariance::analyzeSlopeMethod(allanRes);
const auto& regrRes = AllanVariance::analyzeRegressionMethod(allanRes);

std::cout << AllanVariance::write(slopeRes) << std::endl;
std::cout << AllanVariance::write(regrRes) << std::endl;

## Requires
- Eigen3
- (optional: OpenMP)
