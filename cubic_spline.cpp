#include "cubic_spline.h"
#include <algorithm>

namespace {
    constexpr double kEps = 1e-12;
}

CubicSpline::CubicSpline(const std::vector<Point>& pts) : points_(pts) {
    if (points_.size() < 2) {
        throw std::invalid_argument("Для построения сплайна требуется минимум 2 точки");
    }
    for (size_t i = 1; i < points_.size(); ++i) {
        if (points_[i].x <= points_[i-1].x) {
            throw std::invalid_argument("Узлы должны быть строго упорядочены по возрастанию x");
        }
    }

    const size_t n = points_.size();   // n+1 узлов, n интервалов

    // Шаги
    h_.resize(n - 1);
    for (size_t i = 0; i < n - 1; ++i) {
        h_[i] = points_[i+1].x - points_[i].x;
        if (std::abs(h_[i]) < kEps) {
            throw std::invalid_argument("Шаг между узлами не может быть нулевым");
        }
    }

    // Разделённые разности
    const std::vector<double> delta  = computeDividedDifferences();
    const std::vector<double> delta2 = computeSecondDifferences(delta);
    const std::vector<double> delta3 = computeThirdDifferences(delta2);

    // Трёхдиагональная система для sigma_i, i = 0..n
    std::vector<double> a(n, 0.0), b_coeff(n, 0.0), c_coeff(n, 0.0), d_coeff(n, 0.0);

    // Первое уравнение: -h_0 σ_0 + h_0 σ_1 = h_0² Δ_0^(3)
    a[0]       = 0.0;
    b_coeff[0] = -h_[0];
    c_coeff[0] =  h_[0];
    d_coeff[0] =  h_[0] * h_[0] * delta3[0];

    // Внутренние уравнения
    for (size_t i = 1; i < n - 1; ++i) {
        a[i]       = h_[i-1];
        b_coeff[i] = 2.0 * (h_[i-1] + h_[i]);
        c_coeff[i] = h_[i];
        d_coeff[i] = delta[i] - delta[i-1];
    }

    // Последнее уравнение: h_{n-2} σ_{n-1} - h_{n-2} σ_n = -h_{n-2}² Δ_{n-3}^(3)
    a[n-1]       =  h_[n-2];
    b_coeff[n-1] = -h_[n-2];
    c_coeff[n-1] =  0.0;
    d_coeff[n-1] = -h_[n-2] * h_[n-2] * delta3[n-4];

    sigma_ = solveTridiagonal(a, b_coeff, c_coeff, d_coeff);

    // Коэффициенты на каждом интервале: S(x) = y_i + b_i dx + c_i dx² + d_i dx³
    b_.resize(n - 1);
    c_.resize(n - 1);
    d_.resize(n - 1);
    for (size_t i = 0; i < n - 1; ++i) {
        b_[i] = (points_[i+1].y - points_[i].y) / h_[i] - h_[i] * (sigma_[i+1] + 2.0 * sigma_[i]);
        c_[i] = 3.0 * sigma_[i];
        d_[i] = (sigma_[i+1] - sigma_[i]) / h_[i];
    }
}

std::vector<double> CubicSpline::computeDividedDifferences() const {
    const size_t n = points_.size();
    std::vector<double> delta(n - 1);
    for (size_t i = 0; i < n - 1; ++i) {
        delta[i] = (points_[i+1].y - points_[i].y) / h_[i];
    }
    return delta;
}

std::vector<double> CubicSpline::computeSecondDifferences(const std::vector<double>& delta) const {
    const size_t n = points_.size();
    std::vector<double> delta2(n - 2);
    for (size_t i = 0; i < n - 2; ++i) {
        delta2[i] = (delta[i+1] - delta[i]) / (points_[i+2].x - points_[i].x);
    }
    return delta2;
}

std::vector<double> CubicSpline::computeThirdDifferences(const std::vector<double>& delta2) const {
    const size_t n = points_.size();
    std::vector<double> delta3(n - 3);
    for (size_t i = 0; i < n - 3; ++i) {
        delta3[i] = (delta2[i+1] - delta2[i]) / (points_[i+3].x - points_[i].x);
    }
    return delta3;
}

std::vector<double> CubicSpline::solveTridiagonal(
        const std::vector<double>& a,
        const std::vector<double>& b_coeff,
        const std::vector<double>& c_coeff,
        const std::vector<double>& d_coeff)
{
    const size_t n = a.size();
    if (n == 0) return {};
    if (b_coeff[0] == 0.0) throw std::runtime_error("Прогонка: нулевой диагональный элемент");

    std::vector<double> alpha(n, 0.0), beta(n, 0.0), x(n, 0.0);

    // Прямой ход
    alpha[0] = -c_coeff[0] / b_coeff[0];
    beta[0]  =  d_coeff[0] / b_coeff[0];
    for (size_t i = 1; i < n; ++i) {
        const double denom = b_coeff[i] + a[i] * alpha[i-1];
        if (std::abs(denom) < kEps) {
            throw std::runtime_error("Прогонка: деление на ноль");
        }
        alpha[i] = -c_coeff[i] / denom;
        beta[i]  = (d_coeff[i] - a[i] * beta[i-1]) / denom;
    }

    // Обратный ход
    x[n-1] = beta[n-1];
    for (int i = static_cast<int>(n) - 2; i >= 0; --i) {
        x[i] = alpha[i] * x[i+1] + beta[i];
    }
    return x;
}

size_t CubicSpline::findInterval(double x) const {
    const size_t n = points_.size();
    if (x < points_.front().x - kEps || x > points_.back().x + kEps) {
        throw std::out_of_range("Аргумент находится вне интервала интерполяции");
    }
    // Бинарный поиск интервала
    auto it = std::upper_bound(points_.begin(), points_.end(), x,
        [](double v, const Point& p) { return v < p.x; });
    if (it == points_.begin()) return 0;
    if (it == points_.end()) return n - 2;
    return static_cast<size_t>(std::distance(points_.begin(), it)) - 1;
}

double CubicSpline::interpolate(double x) const {
    const size_t i = findInterval(x);
    const double dx = x - points_[i].x;
    return points_[i].y + b_[i] * dx + c_[i] * dx * dx + d_[i] * dx * dx * dx;
}

double CubicSpline::derivative1(double x) const {
    const size_t i = findInterval(x);
    const double dx = x - points_[i].x;
    return b_[i] + 2.0 * c_[i] * dx + 3.0 * d_[i] * dx * dx;
}

double CubicSpline::derivative2(double x) const {
    const size_t i = findInterval(x);
    const double dx = x - points_[i].x;
    return 2.0 * c_[i] + 6.0 * d_[i] * dx;
}