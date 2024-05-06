#pragma once
#include "rust/cxx.h"
#include <memory>

// Namespace to handle calls into geogram from rust.
namespace GEOGRAM {
    
int16_t det_3d(const ::std::array<double, 3> &a, const ::std::array<double, 3> &b, const ::std::array<double, 3> &c);
int16_t geo_sgn(double x);
int16_t in_circle_2d_SOS(const ::std::array<double, 2> &a, const ::std::array<double, 2> &b, const ::std::array<double, 2> &c, const ::std::array<double, 2> &p);
int16_t in_sphere_3d_SOS(const ::std::array<double, 3> &a, const ::std::array<double, 3> &b, const ::std::array<double, 3> &c, const ::std::array<double, 3> &d, const ::std::array<double, 3> &p);
void initialize();
int16_t orient_2d(const ::std::array<double, 2> &a, const  ::std::array<double, 2> &b, const ::std::array<double, 2> &c);
int16_t orient_2dlifted_SOS(const ::std::array<double, 2> &a, const ::std::array<double, 2> &b, const ::std::array<double, 2> &c, const ::std::array<double, 2> &p, double h_a, double h_b, double h_c, double h_p);
int16_t orient_3d(const ::std::array<double, 3> &a, const  ::std::array<double, 3> &b, const ::std::array<double, 3> &c, const ::std::array<double, 3> &d);
int16_t orient_3dlifted_SOS(const ::std::array<double, 3> &a, const ::std::array<double, 3> &b, const ::std::array<double, 3> &c, const ::std::array<double, 3> &d, const ::std::array<double, 3> &p, double h_a, double h_b, double h_c, double h_d, double h_p);
bool points_are_identical_2d(const ::std::array<double, 2> &p1, const ::std::array<double, 2> &p2);
bool points_are_identical_3d(const ::std::array<double, 3> &p1, const ::std::array<double, 3> &p2);
void show_stats();
void terminate();
} // namespace GEOGRAM