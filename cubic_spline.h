#pragma once

#include <vector>
#include <stdexcept>
#include <cmath>
#include "point.h"

/**
 * @brief Естественный кубический сплайн.
 *
 * Строит кусочно-кубическую интерполирующую функцию S(x) такую, что
 * S(x_i) = y_i, S, S', S'' непрерывны, а S''(x_0) = S''(x_n) = 0.
 * Использует трёхдиагональную систему для параметров sigma_i.
 */
class CubicSpline {
public:
    /**
     * @param pts узлы, строго упорядоченные по возрастанию x
     * @throws std::invalid_argument если точек < 2 или узлы не упорядочены
     */
    explicit CubicSpline(const std::vector<Point>& pts);

    /// Значение сплайна в точке x
    double interpolate(double x) const;

    /// Первая производная S'(x)
    double derivative1(double x) const;

    /// Вторая производная S''(x)
    double derivative2(double x) const;

    /// Коэффициенты сплайна на i-м интервале: S(x)=y_i+b_i dx+c_i dx^2+d_i dx^3
    double b(size_t i) const { return b_[i]; }
    double c(size_t i) const { return c_[i]; }
    double d(size_t i) const { return d_[i]; }
    double sigma(size_t i) const { return sigma_[i]; }

private:
    std::vector<Point>  points_;
    std::vector<double> h_;       ///< шаги h_i = x_{i+1} - x_i
    std::vector<double> sigma_;   ///< параметры sigma_i, i = 0..n
    std::vector<double> b_, c_, d_;

    // Промежуточные вычисления (объявлены private, определены в .cpp)
    std::vector<double> computeDividedDifferences() const;
    std::vector<double> computeSecondDifferences(const std::vector<double>& delta) const;
    std::vector<double> computeThirdDifferences(const std::vector<double>& delta2) const;

    /**
     * @brief Решает трёхдиагональную систему A x = d методом прогонки.
     * a[i]*x[i-1] + b_coeff[i]*x[i] + c_coeff[i]*x[i+1] = d_coeff[i].
     */
    static std::vector<double> solveTridiagonal(
        const std::vector<double>& a,
        const std::vector<double>& b_coeff,
        const std::vector<double>& c_coeff,
        const std::vector<double>& d_coeff);

    /// Индекс интервала [x_i, x_{i+1}], содержащего x (0..n-2)
    size_t findInterval(double x) const;
};