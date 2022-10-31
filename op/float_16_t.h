#ifndef FLOAT_16_T_H
#define FLOAT_16_T_H

#include <iostream>
#include <limits>

constexpr int _v_size = 16; // avx512 / 4bytes = 16 floats
constexpr float infty = std::numeric_limits<float>::infinity();

typedef float float16_t __attribute__ ((vector_size (16 * sizeof(float))));

static float16_t* float16_alloc(std::size_t n) {
    void* tmp = 0;
    if (posix_memalign(&tmp, sizeof(float16_t), sizeof(float16_t) * n)) {
        throw std::bad_alloc();
    }
    return (float16_t*)tmp;
}

constexpr float16_t f16infty {
    infty, infty, infty, infty,
    infty, infty, infty, infty,
    infty, infty, infty, infty,
    infty, infty, infty, infty
};

constexpr float16_t f16zero{
    0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0
};

static inline float hmin16(float16_t vv) {
    float v = infty;
    for (int i = 0; i < 8; ++i) {
        v = std::min(vv[i], v);
    }
    return v;
}

static inline float16_t min8(float16_t x, float16_t y) {
    return x < y ? x : y;
}

#endif