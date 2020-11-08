#include <BicubicSplines.h>
#include <Utility.h>
#include <boost/math/differentiation/finite_difference.hpp>
#include <array>
#include <chrono>
#include <iostream>
#include <random>
#include <Eigen/Dense>
#include <InterpolantBuilder.h>

float func(float x_1, float x_2)
{
    return x_1 * x_1 + x_2 * x_2 + x_1 * x_2;
}

int main(int argc, char *argv[])
{
    constexpr size_t size = 1001;

    auto low = std::array<float, 2>{ -1, -1 };
    auto high = std::array<float, 2>{ 1, 1 };
    auto step_size = std::array<float, 2>{ (high[0] - low[0])/(size-1), (high[0] - low[0])/(size-1) };

    auto def = BicubicSplines::Definition();
    def.f = func;
    def.x_trafo[0] = std::make_unique<LinAxis>(low[0], step_size[0]);
    def.x_trafo[1] = std::make_unique<LinAxis>(low[1], step_size[1]);
    def.nodes = std::array<size_t, 2>{size,size};

    auto builder = InterpolantBuilder<BicubicSplines>(std::move(def));
    auto splines = builder.build();

    auto idef = InterpolationUtility<BicubicSplines>::Definition();
    idef.x_trafo[0] = std::make_unique<LinAxis>(low[0], step_size[0]);
    idef.x_trafo[1] = std::make_unique<LinAxis>(low[1], step_size[1]);

    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<float> dis(-1.0, 1.0);

    /* auto splines = BicubicSplines(y, dydx1, dydx2, d2ydx1dx2); */
    auto inter = InterpolationUtility<BicubicSplines>(std::move(splines), std::move(idef));

    auto point = std::array<float, 2>();
    auto res = 0.f;
    auto start = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < 1'000'000; ++i) {
        point[0] = dis(gen);
        point[1] = dis(gen);

        res = inter.evaluate(point);
    }
    auto stop = std::chrono::high_resolution_clock::now();

    std::cout << point[0] << " * " << point[0] << " + "<< point[1] << " * " << point[1] << " + " << point[0] << " * "<< point[1]  << " = " << res << std::endl;

    std::cout << "cubic splien interpolation takes "
        << std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count()
        << " micro seconds"
        << std::endl;


    return 0;
}
