#include "include/geogram_predicates_psm/Predicates_psm.h"

#include <array>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {

enum Predicate : std::uint8_t {
    ORIENT_2D = 0,
    ORIENT_3D = 1,
    DET_3D = 2,
    DET_4D = 3,
    DET_COMPARE_4D = 4,
    DOT_3D = 5,
    DOT_COMPARE_3D = 6,
    ALIGNED_3D = 7,
    IDENTICAL_2D = 8,
    IDENTICAL_3D = 9,
    COLINEAR_3D = 10,
    IN_CIRCLE_2D_SOS = 11,
    IN_CIRCLE_3D_SOS = 12,
    IN_SPHERE_3D_SOS = 13,
    SIDE1_SOS = 14,
    SIDE2_SOS = 15,
    SIDE3_SOS = 16,
    SIDE4_SOS = 17,
    SIDE4_3D = 18,
    SIDE4_3D_SOS = 19,
    SIDE3_3D_LIFTED = 20,
    SIDE3_3D_LIFTED_SOS = 21,
    IN_CIRCLE_3D_LIFTED = 22,
    IN_CIRCLE_3D_LIFTED_SOS = 23,
    ORIENT_2D_LIFTED_SOS = 24,
    ORIENT_3D_LIFTED = 25,
    ORIENT_3D_LIFTED_SOS = 26,
};

std::uint8_t read_byte(bool allow_eof = false) {
    const int value = std::cin.get();
    if(value == std::char_traits<char>::eof()) {
        if(allow_eof) {
            return 0xff;
        }
        throw std::runtime_error("truncated differential-test input");
    }
    return static_cast<std::uint8_t>(value);
}

double read_f64() {
    std::uint64_t bits = 0;
    for(unsigned int shift = 0; shift != 64; shift += 8) {
        bits |= static_cast<std::uint64_t>(read_byte()) << shift;
    }
    double value;
    static_assert(sizeof(value) == sizeof(bits), "binary64 is required");
    std::memcpy(&value, &bits, sizeof(value));
    return value;
}

int evaluate_case(std::uint8_t predicate) {
    const std::uint8_t dim = read_byte();
    const std::uint8_t point_count = read_byte();
    const std::uint8_t query_count = read_byte();
    const std::uint8_t height_count = read_byte();
    if(dim > 8 || point_count > 9 || query_count > 4 || height_count > 5) {
        throw std::runtime_error("invalid differential-test record shape");
    }

    std::array<std::uint8_t, 9> keys{};
    std::array<bool, 9> key_seen{};
    for(std::uint8_t i = 0; i < point_count; ++i) {
        keys[i] = read_byte();
        if(keys[i] >= point_count || key_seen[keys[i]]) {
            throw std::runtime_error("SOS keys are not a permutation");
        }
        key_seen[keys[i]] = true;
    }

    // All point pointers refer into this single array. Pointer order is therefore
    // defined by the key-sized stride and matches Rust's stable key order.
    std::array<double, 9 * 8> point_storage{};
    for(std::uint8_t i = 0; i < point_count; ++i) {
        for(std::uint8_t coordinate = 0; coordinate < dim; ++coordinate) {
            point_storage[static_cast<std::size_t>(keys[i]) * 8 + coordinate] = read_f64();
        }
    }
    std::array<const double*, 9> p{};
    for(std::uint8_t i = 0; i < point_count; ++i) {
        p[i] = point_storage.data() + static_cast<std::size_t>(keys[i]) * 8;
    }

    std::array<double, 4 * 8> query_storage{};
    std::array<const double*, 4> q{};
    for(std::uint8_t i = 0; i < query_count; ++i) {
        q[i] = query_storage.data() + static_cast<std::size_t>(i) * 8;
        for(std::uint8_t coordinate = 0; coordinate < dim; ++coordinate) {
            query_storage[static_cast<std::size_t>(i) * 8 + coordinate] = read_f64();
        }
    }

    std::array<double, 5> h{};
    for(std::uint8_t i = 0; i < height_count; ++i) {
        h[i] = read_f64();
    }

    using namespace GEO::PCK;
    switch(static_cast<Predicate>(predicate)) {
    case ORIENT_2D:
        return orient_2d(p[0], p[1], p[2]);
    case ORIENT_3D:
        return orient_3d(p[0], p[1], p[2], p[3]);
    case DET_3D:
        return det_3d(p[0], p[1], p[2]);
    case DET_4D:
        return det_4d(p[0], p[1], p[2], p[3]);
    case DET_COMPARE_4D:
        return det_compare_4d(p[0], p[1], p[2], p[3], p[4]);
    case DOT_3D:
        return dot_3d(p[0], p[1], p[2]);
    case DOT_COMPARE_3D:
        return dot_compare_3d(p[0], p[1], p[2]);
    case ALIGNED_3D:
        return aligned_3d(p[0], p[1], p[2]) ? 1 : 0;
    case IDENTICAL_2D:
        return points_are_identical_2d(p[0], p[1]) ? 1 : 0;
    case IDENTICAL_3D:
        return points_are_identical_3d(p[0], p[1]) ? 1 : 0;
    case COLINEAR_3D:
        return points_are_colinear_3d(p[0], p[1], p[2]) ? 1 : 0;
    case IN_CIRCLE_2D_SOS:
        return in_circle_2d_SOS(p[0], p[1], p[2], p[3]);
    case IN_CIRCLE_3D_SOS:
        return in_circle_3d_SOS(p[0], p[1], p[2], p[3]);
    case IN_SPHERE_3D_SOS:
        return in_sphere_3d_SOS(p[0], p[1], p[2], p[3], p[4]);
    case SIDE1_SOS:
        return side1_SOS(p[0], p[1], q[0], dim);
    case SIDE2_SOS:
        return side2_SOS(p[0], p[1], p[2], q[0], q[1], dim);
    case SIDE3_SOS:
        return side3_SOS(p[0], p[1], p[2], p[3], q[0], q[1], q[2], dim);
    case SIDE4_SOS:
        return side4_SOS(p[0], p[1], p[2], p[3], p[4], q[0], q[1], q[2], q[3], dim);
    case SIDE4_3D:
        return side4_3d(p[0], p[1], p[2], p[3], p[4]);
    case SIDE4_3D_SOS:
        return side4_3d_SOS(p[0], p[1], p[2], p[3], p[4]);
    case SIDE3_3D_LIFTED:
        return side3_3dlifted_SOS(
            p[0], p[1], p[2], p[3], h[0], h[1], h[2], h[3],
            q[0], q[1], q[2], false
        );
    case SIDE3_3D_LIFTED_SOS:
        return side3_3dlifted_SOS(
            p[0], p[1], p[2], p[3], h[0], h[1], h[2], h[3],
            q[0], q[1], q[2], true
        );
    case IN_CIRCLE_3D_LIFTED:
        return in_circle_3dlifted_SOS(
            p[0], p[1], p[2], p[3], h[0], h[1], h[2], h[3], false
        );
    case IN_CIRCLE_3D_LIFTED_SOS:
        return in_circle_3dlifted_SOS(
            p[0], p[1], p[2], p[3], h[0], h[1], h[2], h[3], true
        );
    case ORIENT_2D_LIFTED_SOS:
        return orient_2dlifted_SOS(
            p[0], p[1], p[2], p[3], h[0], h[1], h[2], h[3]
        );
    case ORIENT_3D_LIFTED:
        return orient_3dlifted(
            p[0], p[1], p[2], p[3], p[4], h[0], h[1], h[2], h[3], h[4]
        );
    case ORIENT_3D_LIFTED_SOS:
        return orient_3dlifted_SOS(
            p[0], p[1], p[2], p[3], p[4], h[0], h[1], h[2], h[3], h[4]
        );
    }
    throw std::runtime_error("unknown differential-test predicate");
}

} // namespace

int main() {
    try {
        GEO::PCK::initialize();
        GEO::PCK::set_SOS_mode(GEO::PCK::SOS_ADDRESS);
        while(true) {
            const std::uint8_t predicate = read_byte(true);
            if(predicate == 0xff) {
                break;
            }
            const int result = evaluate_case(predicate);
            if(result < -1 || result > 1) {
                throw std::runtime_error("oracle returned an invalid result");
            }
            std::cout.put(static_cast<char>(static_cast<std::int8_t>(result)));
        }
        GEO::PCK::terminate();
        return std::cout ? 0 : 3;
    } catch(const std::exception& error) {
        std::cerr << "C++ differential oracle: " << error.what() << '\n';
        return 2;
    }
}
