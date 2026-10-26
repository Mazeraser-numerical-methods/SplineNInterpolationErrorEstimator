#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <string>

#include "point.h"
#include "cubic_spline.h"
#include "interpolation_error_estimator.h"

static void printSplineTable(const CubicSpline& spline,
                             double x0, double x1, double step)
{
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "     x          S(x)          S'(x)         S''(x)\n";
    std::cout << "--------------------------------------------------\n";
    for (double x = x0; x <= x1 + 1e-9; x += step) {
        std::cout << std::setw(8)  << x
                  << std::setw(14) << spline.interpolate(x)
                  << std::setw(14) << spline.derivative1(x)
                  << std::setw(14) << spline.derivative2(x)
                  << '\n';
    }
    std::cout << '\n';
}

int main() {
    try {
        // ---- Пример 6 из лекции 4 ----
        std::vector<Point> pts = {
            {1.0, 2.0}, {3.0, 4.0}, {6.0, 7.0},
            {7.0, 7.0}, {8.0, 6.0}, {9.0, 5.0}
        };

        std::cout << "=== Сплайн-интерполяция ===\n";
        CubicSpline spline(pts);

        std::cout << "Значения сплайна на [1, 9]:\n";
        printSplineTable(spline, 1.0, 9.0, 0.5);

        std::cout << "Проверка прохождения через узлы:\n";
        for (const auto& p : pts) {
            std::cout << "  S(" << p.x << ") = " << spline.interpolate(p.x)
                      << "  (ожидается " << p.y << ")\n";
        }
        std::cout << '\n';

        // ---- Оценка погрешности ----
        std::cout << "=== Оценка погрешности интерполяционных формул ===\n";

        {
            double q = 0.2;
            double d4 = 1e-4;
            int n = 3;
            double e = InterpolationErrorEstimator::newtonForwardError(q, d4, n);
            std::cout << "Ньютон (вперёд): q=" << q << ", n=" << n
                      << ", |R| ~ " << e << '\n';
        }
        {
            double q = -0.7;
            double d4 = 1e-4;
            int n = 3;
            double e = InterpolationErrorEstimator::newtonBackwardError(q, d4, n);
            std::cout << "Ньютон (назад):  q=" << q << ", n=" << n
                      << ", |R| ~ " << e << '\n';
        }
        {
            double q = 0.1;
            double d3 = 1e-3;
            int n = 1;
            double e = InterpolationErrorEstimator::stirlingError(q, d3, n);
            std::cout << "Стирлинг:        q=" << q << ", n=" << n
                      << ", |R| ~ " << e << '\n';
        }
        {
            double q = 0.5;
            double d4 = 1e-4;
            int n = 2;
            double e = InterpolationErrorEstimator::besselError(q, d4, n);
            std::cout << "Бессель:         q=" << q << ", n=" << n
                      << ", |R| ~ " << e << '\n';
        }
        std::cout << '\n';

        // ---- Пример 7: монотонные данные, сплайн vs полином ----
        std::cout << "=== Пример 7: монотонно убывающие данные ===\n";
        std::vector<Point> bad = {
            {4.0, 3.73}, {4.5, 2.37}, {5.5, 2.14},
            {6.0, 0.55}, {7.0, 0.09}, {7.5, 0.03}
        };
        CubicSpline monotonic(bad);
        std::cout << "Сплайн сохраняет монотонность:\n";
        std::cout << "     x          S(x)\n";
        for (double x = 4.0; x <= 7.5 + 1e-9; x += 0.1) {
            std::cout << std::setw(8) << x
                      << std::setw(14) << monotonic.interpolate(x) << '\n';
        }
        std::cout << '\n';

        // ---- Сплайн для синуса ----
        std::cout << "=== Сплайн-интерполяция синуса ===\n";
        std::vector<Point> sinPts;
        for (int i = 0; i < 8; ++i) {
            double x = i * 0.5;
            sinPts.push_back({x, std::sin(x)});
        }
        CubicSpline sinSpline(sinPts);
        std::cout << "     x         S(x)          sin(x)        погрешность\n";
        for (double x = 0.0; x <= 3.5 + 1e-9; x += 0.1) {
            double s = sinSpline.interpolate(x);
            double t = std::sin(x);
            std::cout << std::setw(8)  << x
                      << std::setw(14) << s
                      << std::setw(14) << t
                      << std::setw(14) << std::abs(s - t) << '\n';
        }

    } catch (const std::exception& e) {
        std::cerr << "Ошибка: " << e.what() << '\n';
        return 1;
    }
    return 0;
}