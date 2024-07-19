#include "geogram_predicates/include/geogram_ffi.h"
#include "geogram_predicates/include/geogram_predicates_psm/Predicates_psm.h"
#include <algorithm>
#include <functional>
#include <iomanip>
#include <set>
#include <string>
#include <unordered_map>

namespace GEOGRAM {

int16_t det_3d(const ::std::array<double, 3> &a, const ::std::array<double, 3> &b, const ::std::array<double, 3> &c) {
    GEO::Sign det = GEO::PCK::det_3d(a.data(), b.data(), c.data());
    return det;
}

int16_t det_4d(const ::std::array<double, 4> &a, const ::std::array<double, 4> &b, const ::std::array<double, 4> &c, const ::std::array<double, 4> &d) {
    GEO::Sign det = GEO::PCK::det_4d(a.data(), b.data(), c.data(), d.data());
    return det;
}

int16_t dot_3d(const ::std::array<double, 3> &a, const ::std::array<double, 3> &b, const ::std::array<double, 3> &c) {
    GEO::Sign dot = GEO::PCK::dot_3d(a.data(), b.data(), c.data());
    return dot;
}

int16_t geo_sgn(double x) {
    return GEO::geo_sgn(x);
};

int16_t in_circle_2d_SOS(const ::std::array<double, 2> &a, const ::std::array<double, 2> &b, const ::std::array<double, 2> &c, const ::std::array<double, 2> &p) {
    GEO::Sign is_in = GEO::PCK::in_circle_2d_SOS(a.data(), b.data(), c.data(), p.data());
    return is_in;
}

int16_t in_sphere_3d_SOS(const ::std::array<double, 3> &a, const ::std::array<double, 3> &b, const ::std::array<double, 3> &c, const ::std::array<double, 3> &d, const ::std::array<double, 3> &p) {
    GEO::Sign is_in = GEO::PCK::in_sphere_3d_SOS(a.data(), b.data(), c.data(), d.data(), p.data());
    return is_in;
}

void initialize() {
    GEO::PCK::initialize();
}

int16_t orient_2d(const ::std::array<double, 2> &a, const ::std::array<double, 2> &b, const ::std::array<double, 2> &c) {
    GEO::Sign orientation = GEO::PCK::orient_2d(a.data(), b.data(), c.data());
    return orientation;
}

int16_t orient_2dlifted_SOS(const ::std::array<double, 2> &a, const ::std::array<double, 2> &b, const ::std::array<double, 2> &c, const ::std::array<double, 2> &p, double h_a, double h_b, double h_c, double h_p) {
    GEO::Sign regularity = GEO::PCK::orient_2dlifted_SOS(a.data(), b.data(), c.data(), p.data(), h_a, h_b, h_c, h_p);
    return regularity;
}

int16_t orient_3d(const ::std::array<double, 3> &a, const  ::std::array<double, 3> &b, const ::std::array<double, 3> &c, const ::std::array<double, 3> &d) {
    GEO::Sign orientation = GEO::PCK::orient_3d(a.data(), b.data(), c.data(), d.data());
    return orientation;
}

int16_t orient_3d_inexact(const ::std::array<double, 3> &a, const  ::std::array<double, 3> &b, const ::std::array<double, 3> &c, const ::std::array<double, 3> &d) {
    GEO::Sign orientation = GEO::PCK::orient_3d_inexact(a.data(), b.data(), c.data(), d.data());
    return orientation;
}

int16_t orient_3dlifted_SOS(const ::std::array<double, 3> &a, const ::std::array<double, 3> &b, const ::std::array<double, 3> &c, const ::std::array<double, 3> &d, const ::std::array<double, 3> &p, double h_a, double h_b, double h_c, double h_d, double h_p) {
    GEO::Sign regularity = GEO::PCK::orient_3dlifted_SOS(a.data(), b.data(), c.data(), d.data(), p.data(), h_a, h_b, h_c, h_d, h_p);
    return regularity;
}

bool points_are_identical_2d(const ::std::array<double, 2> &p1, const ::std::array<double, 2> &p2) {
    return GEO::PCK::points_are_identical_2d(p1.data(), p2.data());
}

bool points_are_identical_3d(const ::std::array<double, 3> &p1, const ::std::array<double, 3> &p2) {
    return GEO::PCK::points_are_identical_3d(p1.data(), p2.data());
}

void show_stats() {
    GEO::PCK::show_stats();
}

void terminate() {
    GEO::PCK::terminate();
}

} // namespace GEOGRAM