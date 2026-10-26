#include "interpolation_error_estimator.h"
#include <stdexcept>

namespace {
    // Факториал через tgamma: n! = Γ(n+1)
    inline double factorial(int n) { return std::tgamma(n + 1); }
}

double InterpolationErrorEstimator::newtonForwardError(double q, double delta_n_plus_1, int n) {
    if (n < 0) throw std::invalid_argument("Степень полинома должна быть >= 0");
    double product = 1.0;
    for (int i = 0; i <= n; ++i) {
        product *= (q - i);            // q(q-1)...(q-n)
    }
    return std::abs(product * delta_n_plus_1 / factorial(n + 1));
}

double InterpolationErrorEstimator::newtonBackwardError(double q, double delta_n_plus_1, int n) {
    if (n < 0) throw std::invalid_argument("Степень полинома должна быть >= 0");
    double product = 1.0;
    for (int i = 0; i <= n; ++i) {
        product *= (q + i);            // q(q+1)...(q+n)
    }
    return std::abs(product * delta_n_plus_1 / factorial(n + 1));
}

double InterpolationErrorEstimator::stirlingError(double q, double delta_odd, int n) {
    if (n < 0) throw std::invalid_argument("Степень полинома должна быть >= 0");
    // q(q²-1)(q²-4)...(q²-n²)
    double product = q;
    for (int i = 1; i <= n; ++i) {
        product *= (q * q - static_cast<double>(i) * i);
    }
    return std::abs(product * delta_odd / (2.0 * factorial(2 * n + 1)));
}

double InterpolationErrorEstimator::besselError(double q, double delta_even, int n) {
    if (n < 0) throw std::invalid_argument("Степень полинома должна быть >= 0");
    // q(q²-1)...(q²-n²)(q-n-1) = q(q-1)(q+1)...(q-n)(q+n)  (эквивалентно)
    double product = q;
    for (int i = 1; i <= n; ++i) {
        product *= (q * q - static_cast<double>(i) * i);
    }
    product *= (q - n - 1);
    return std::abs(product * delta_even / (2.0 * factorial(2 * n + 2)));
}