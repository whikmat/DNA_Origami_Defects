/**
 * @brief SIMD functions and Boost.SIMD wrapper for dynamic programs
 *
 * @file SIMD.h
 * @author Mark Fornace
 * @date 2018-05-31
 */
#pragma once
#include "../iteration/Range.h"
#include "../common/Error.h"
#include "../reflect/Print.h"

#include <boost/align/aligned_allocator.hpp>

/******************************************************************************************/

#ifndef NUPACK_NO_SIMD
#    include <simdpp/simd.h>
#endif

// If the eval is not used, it might be faster on some FMA architectures.
// But it gets unsafe due to returning internal references so the file with LIBSIMDPP_SIMD_EXPR_H 
// defined has to be changed so all the members are not references.
// There isn't any real difference I can see on arm64 mac so I opted against modifying it.
// #define NUPACK_SIMD_EVAL(x) x
#define NUPACK_SIMD_EVAL(x) (x).eval()

/******************************************************************************************/

namespace nupack::math {

    template <class To, class From, NUPACK_IF(sizeof(To) == sizeof(From))>
    To bit_cast(From const &src) noexcept {
        To dst;
        std::memcpy(&dst, &src, sizeof(To));
        return dst;
    }

    static_assert(std::numeric_limits<float>::max_exponent == 128);
    static_assert(std::numeric_limits<double>::max_exponent == 1024);

    // ldexp: build float of 2^e, then multiply by t. Denormals are not returned.
    // the min max is necessary so that overflow is captured as inf and underflow as 0
    // Incurred operations: 1 each of *, max, min, +, <<

    template <class T, NUPACK_IF(is_same<T, float>)>
    float ldexp(T t, std::int32_t e) noexcept {
        using I = std::int32_t;
        return t * bit_cast<float>(std::max<I>(I(0), std::min<I>(255, e + 127)) << 23);
    }

    template <class T, NUPACK_IF(is_same<T, double>)>
    inline double ldexp(T t, std::int64_t e) noexcept {
        using I = std::int64_t;
        return t * bit_cast<double>(std::max<I>(I(0), std::min<I>(2047, e + 1023)) << 52);
    }    

#ifndef NUPACK_NO_SIMD
    template <class T>
    using scalar_type = std::decay_t<decltype(simdpp::reduce_add(std::declval<T>()))>;

    template <class T, class E, NUPACK_IF(is_same<scalar_type<T>, float>)>
    inline auto ldexp(T const &t, E const &e) noexcept {
        using F = decltype(simdpp::to_float32(t).eval());
        using I = std::int32_t;
        return NUPACK_SIMD_EVAL(t * simdpp::bit_cast<F>(simdpp::shift_l<23>(
            simdpp::max(I(0), simdpp::min(I(255), e + I(127))))));
    }

    template <class T, class E, NUPACK_IF(is_same<scalar_type<T>, double>)>
    inline auto ldexp(T const &t, E const &e) noexcept {
        using F = decltype(simdpp::to_float64(t).eval());
        using I = std::int64_t;
        return NUPACK_SIMD_EVAL(t * simdpp::bit_cast<F>(simdpp::shift_l<52>(
            simdpp::max(I(0), simdpp::min(I(2047), e + I(1023))))));
    }
#endif
}

/******************************************************************************************/

namespace nupack::simd {

/******************************************************************************************/

/* Chunk is used as an element accessor, i.e. array[Chunk()]. It tells the container
   to return N elements starting at a given position  */
template <int N>
struct Chunk {
    static constexpr auto length = N;
    int value;
    explicit constexpr Chunk(int v) : value(v) {}

    friend auto operator-(Chunk c, int i) noexcept {return Chunk(c.value-i);}
    friend auto operator-(int i, Chunk c) noexcept {return Chunk<-N>(i-c.value);}
    friend auto operator+(int i, Chunk c) noexcept {return Chunk(c.value+i);}
    friend auto operator+(Chunk c, int i) noexcept {return Chunk(c.value+i);}
};

NUPACK_DEFINE_TEMPLATE(is_chunk, Chunk, int);

/******************************************************************************************/

template <class T, class SFINAE=void>
struct SingleDispatch;

#define NUPACK_TMP(name) \
    struct name##_t { \
        template <class T> \
        auto operator()(T &&t) const noexcept { \
            return SingleDispatch<no_qual<T>>::name(std::forward<T>(t)); \
        } \
    }; \
    static constexpr name##_t name{};

NUPACK_TMP(ifrexp);
NUPACK_TMP(reciprocal);
NUPACK_TMP(negate);
NUPACK_TMP(reduce_sum);
NUPACK_TMP(reduce_min);

#undef NUPACK_TMP

template <class T>
struct SingleDispatch<T, std::enable_if_t<std::is_arithmetic_v<T>>> {
    static T reciprocal(T const &t) noexcept {return T(1) / t;}
    static T negate(T const &t) noexcept {return -t;}
    static auto ifrexp(T const &t) noexcept {
        int e;
        T m = std::frexp(t, &e);
        return std::pair<T, int_of_size<sizeof(T)>>(m, e);
    }
    static constexpr T const &reduce_sum(T const &t) noexcept {return t;}
    static constexpr T const &reduce_min(T const &t) noexcept {return t;}
};

/******************************************************************************************/

template <class T, class U, class SFINAE=void>
struct DoubleDispatch;

#define NUPACK_TMP(name) \
    struct name##_t { \
        template <class T, class U> \
        auto operator()(T &&t, U &&u) const noexcept { \
            return DoubleDispatch<no_qual<T>, no_qual<U>>::name(std::forward<T>(t), std::forward<U>(u)); \
        } \
    }; \
    static constexpr name##_t name{};

NUPACK_TMP(load);
NUPACK_TMP(store);
NUPACK_TMP(plus);
NUPACK_TMP(min);
NUPACK_TMP(max);
NUPACK_TMP(multiplies);
NUPACK_TMP(ldexp);

#undef NUPACK_TMP

/******************************************************************************************/

template <class T>
struct DoubleDispatch<T, T, std::enable_if_t<std::is_arithmetic_v<T>>> {
    static T plus(T const &t, T const &u) noexcept {return t + u;}
    static T multiplies(T const &t, T const &u) noexcept {return t * u;}
    static T min(T const &t, T const &u) noexcept {return std::min(t, u);}
    static T max(T const &t, T const &u) noexcept {return std::max(t, u);}
};

template <class T>
struct DoubleDispatch<T, int_of_size<sizeof(T)>, std::enable_if_t<std::is_floating_point_v<T>>> {
    static T ldexp(T const &t, int_of_size<sizeof(T)> e) noexcept {return math::ldexp(t, e);}
};

/******************************************************************************************/

#ifndef NUPACK_NO_SIMD
template <class T, int N>
using simd_type = if_t<is_same<T, float>, simdpp::float32<N>,
                    if_t<is_same<T, double>, simdpp::float64<N>,
                      if_t<is_same<T, std::int64_t>, simdpp::int64<N>,
                        if_t<is_same<T, std::int32_t>, simdpp::int32<N>, void>>>>;

NUPACK_DETECT(is_simd, decltype(std::declval<T>().eval()));

/******************************************************************************************/

template <class T, std::size_t ...Is>
std::array<math::scalar_type<T>, sizeof...(Is)> array_impl(T const &t, std::index_sequence<Is...>) {return {simdpp::extract<Is>(t)...};}

template <class T>
auto array_impl(T const &t) {return array_impl(t, std::make_index_sequence<T::length>());}
template <class T>
auto array(T const &t) {return array_impl(t.eval());}

/******************************************************************************************/

template <class T>
struct SingleDispatch<T, std::enable_if_t<traits::is_simd<T>>> {
    static auto reduce_sum(T const &t) noexcept {return simdpp::reduce_add(t);}
    static auto reduce_min(T const &t) noexcept {return simdpp::reduce_min(t);}
    static auto negate(T const &t) noexcept {return NUPACK_SIMD_EVAL(simdpp::neg(t));}
};

/******************************************************************************************/

template <class T, int N>
struct DoubleDispatch<Chunk<N>, T *, std::enable_if_t<(N >= 1)>> {
    static auto load(Ignore, T *t) noexcept {return simd_type<std::remove_cv_t<T>, N>(simdpp::load_u(t));}
};

// Loading reversed chunk includes the pointer t but goes backwards from there
// Should experiment to see if the following permute options are actually faster.

// template <class T>
// struct DoubleDispatch<Chunk<-4>, T *> {
//     using O = simd_type<std::remove_cv_t<T>, 4>;
//     static O load(Ignore, T *t) noexcept {
//         return simdpp::permute4<3,2,1,0>(O(simdpp::load_u(t+(1-4))));
//     }
// };

// template <class T>
// struct DoubleDispatch<Chunk<-2>, T *> {
//     using O = simd_type<std::remove_cv_t<T>, 2>;

//     static O load(Ignore, T *t) noexcept {return simdpp::permute2<1,0>(O(simdpp::load_u(t+(1-2))));}
// };

template <class T>
struct DoubleDispatch<Chunk<-1>, T *> {
    using O = simd_type<std::remove_cv_t<T>, 1>;

    static O load(Ignore, T *t) noexcept {return O(simdpp::load_u(t));}
};


template <int N, class T>
struct DoubleDispatch<Chunk<N>, T *, std::enable_if_t<N <= -1>> {
    using O = simd_type<std::remove_cv_t<T>, -N>;

    template <std::size_t ...Is>
    static O impl(T *t, std::index_sequence<Is...>) noexcept {
        if constexpr(std::is_integral_v<T> && std::is_signed_v<T>) {
            return simdpp::make_int(t[-int(Is)]...);
        } else if constexpr(std::is_integral_v<T> && !std::is_signed_v<T>) {
            return simdpp::make_uint(t[-int(Is)]...);
        } else if constexpr(std::is_floating_point_v<T>) {
            return simdpp::make_float(t[-int(Is)]...);
        } else {
            static_assert(T::no_load_possible);
        }
        
    }

    static O load(Ignore, T *t) noexcept {return impl(t, std::make_index_sequence<-N>());}
};

template <class X, class T>
struct DoubleDispatch<X, T *, std::enable_if_t<traits::is_simd<X>>> {
    static void store(X const &x, T *t) noexcept {
        simdpp::store(t, x);
    }
};

/******************************************************************************************/

template <class T, class U>
struct DoubleDispatch<T, U, std::enable_if_t<traits::is_simd<T> && traits::is_simd<U>>> {
    static auto plus(T const &t, U const &u) noexcept {return NUPACK_SIMD_EVAL(simdpp::add(t, u));}
    static auto multiplies(T const &t, U const &u) noexcept {return NUPACK_SIMD_EVAL(simdpp::mul(t, u));}
    static auto min(T const &t, U const &u) noexcept {return NUPACK_SIMD_EVAL(simdpp::min(t, u));}
    static auto max(T const &t, U const &u) noexcept {return NUPACK_SIMD_EVAL(simdpp::max(t, u));}
    static auto ldexp(T const &t, U const &e) noexcept {return math::ldexp(t, e);}
};

/******************************************************************************************/

template <class T, class U>
struct DoubleDispatch<T, U, std::enable_if_t<traits::is_simd<T> && std::is_arithmetic_v<U>>> {
    static auto plus(T const &t, U const &u) noexcept {return NUPACK_SIMD_EVAL(simdpp::add(t, u));}
    static auto multiplies(T const &t, U const &u) noexcept {return NUPACK_SIMD_EVAL(simdpp::mul(t, u));}
    static auto min(T const &t, U const &u) noexcept {return NUPACK_SIMD_EVAL(simdpp::min(t, u));}
};

template <class T, class U>
struct DoubleDispatch<T, U, std::enable_if_t<traits::is_simd<U> && std::is_arithmetic_v<T>>> {
    static auto plus(T const &t, U const &u) noexcept {return NUPACK_SIMD_EVAL(simdpp::add(t, u));}
    static auto multiplies(T const &t, U const &u) noexcept {return NUPACK_SIMD_EVAL(simdpp::mul(t, u));}
    static auto min(T const &t, U const &u) noexcept {return NUPACK_SIMD_EVAL(simdpp::min(t, u));}
};

template <class T> struct OptimalSize {};
template <> struct OptimalSize<float> {static constexpr auto value = SIMDPP_FAST_FLOAT32_SIZE;};
template <> struct OptimalSize<double> {static constexpr auto value = SIMDPP_FAST_FLOAT64_SIZE;};
template <> struct OptimalSize<std::int32_t> {static constexpr auto value = SIMDPP_FAST_FLOAT32_SIZE;};
template <> struct OptimalSize<std::int64_t> {static constexpr auto value = SIMDPP_FAST_FLOAT64_SIZE;};

template <class T>
static constexpr auto optimal_size = OptimalSize<T>::value;

template <class T>
static constexpr std::uintptr_t alignment = sizeof(T) * optimal_size<T>;

template <class T>
using allocator = boost::alignment::aligned_allocator<T, alignment<T>>;

#else

template <class T>
static constexpr std::uintptr_t alignment = alignof(T);

template <class T>
using allocator = std::allocator<T>;

#endif

/******************************************************************************************/

// Compile-time zero object. Avoids computation when an operand is known to be 0.
struct Zero {
    template <class T>
    friend T const &operator+(T const &t, Zero) {return t;}

    template <class T>
    friend T const &operator+(Zero, T const &t) {return t;}

    friend Zero operator+(Zero, Zero) {return {};}

    friend False operator<(Zero, Zero) {return {};}

    template <class T>
    friend Zero operator*(T const &t, Zero) {return {};}

    template <class T>
    friend Zero operator*(Zero, T const &t) {return {};}

    friend Zero operator*(Zero, Zero) {return {};}

    template <class T>
    friend T const &operator-(T const &t, Zero) {return t;}

    template <class T>
    friend T operator-(Zero, T const &t) {return -t;}

    friend Zero operator-(Zero, Zero) {return {};}

    Zero operator-() const {return {};}

    constexpr operator std::int32_t() const {return 0;}
    constexpr operator std::int64_t() const {return 0;}
};

NUPACK_UNARY_FUNCTOR(always_zero, Zero());

/******************************************************************************************/

template <class T>
struct DoubleDispatch<T, Zero> {
    static T const &ldexp(T const &t, Zero) noexcept {return t;}
    static T const &plus(T const &t, Zero) noexcept {return t;}
};

template <>
struct DoubleDispatch<Zero, Zero> {
    static constexpr Zero plus(Zero, Zero) noexcept {return {};}
    static constexpr Zero min(Zero, Zero) noexcept {return {};}
    static constexpr Zero max(Zero, Zero) noexcept {return {};}
};

template <>
struct SingleDispatch<Zero> {
    static constexpr Zero negate(Zero) noexcept {return {};}
};

/******************************************************************************************/

/// Perform a map-reduce operation with SIMD, where the output is modified in place
template <class R, class D, class M, class Op>
auto map_reduce(R reduce, D domain, M map, Op op) noexcept {
    auto it = begin_of(domain);
    auto ret = map(*it++);

#ifndef NUPACK_NO_SIMD
    constexpr auto Z = optimal_size<decltype(ret)>;
    if (it + Z <= end_of(domain)) {
        auto sum = map(Chunk<Z>(*it)).eval();
        for (it += Z; it + Z <= end_of(domain); it += Z) sum = reduce(sum, map(Chunk<Z>(*it)));
        ret = reduce(ret, op(sum));
    }
#endif

    for (; it != end_of(domain); ++it) ret = reduce(ret, map(*it));
    return ret;
}

/******************************************************************************************/

/// Variadic max function
template <class T, class U, class ...Ts>
constexpr auto max_of(T const &t, U const &u, Ts const &...ts) noexcept {
    if constexpr(sizeof...(Ts) == 0) {
        return max(t, u);
    } else {
        return max_of(max(t, u), ts...);
    }
}

/******************************************************************************************/

}
