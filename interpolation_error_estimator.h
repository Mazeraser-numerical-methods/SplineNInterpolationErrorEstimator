#pragma once

#include <cmath>

/**
 * @brief Оценки остаточного члена интерполяционных формул (лекция 5).
 *
 * Все оценки используют приближение производной (n+1)-го порядка
 * через конечные разности: M_{n+1} ≈ Δ^{n+1} y / h^{n+1}.
 */
class InterpolationErrorEstimator {
public:
    /**
     * @brief |R_1| ≈ |q(q-1)...(q-n) / (n+1)! · Δ^{n+1} y_0|
     * @param q = (x - x0)/h
     * @param delta_n_plus_1 конечная разность (n+1)-го порядка в точке x0
     * @param n степень полинома (n >= 0)
     */
    static double newtonForwardError(double q, double delta_n_plus_1, int n);

    /**
     * @brief |R_2| ≈ |q(q+1)...(q+n) / (n+1)! · Δ^{n+1} y_n|
     * @param q = (x - xn)/h
     */
    static double newtonBackwardError(double q, double delta_n_plus_1, int n);

    /**
     * @brief Оценка для формулы Стирлинга:
     * |R_S| ≈ |q(q²-1)(q²-4)...(q²-n²) / (2·(2n+1)!) · Δ^{2n+1}|
     */
    static double stirlingError(double q, double delta_odd, int n);

    /**
     * @brief Оценка для формулы Бесселя:
     * |R_B| ≈ |q(q²-1)...(q²-n²)(q-n-1) / (2·(2n+2)!) · Δ^{2n+2}|
     */
    static double besselError(double q, double delta_even, int n);
};