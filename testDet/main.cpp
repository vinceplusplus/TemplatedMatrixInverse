//
//  main.cpp
//  testDet
//
//  Created by Vincent Cheung on 6/4/16.
//  Copyright © 2016 test. All rights reserved.
//

#include <iostream>
#include <vector>
#include <cmath>
#include <type_traits>

using namespace std;

// TODO: consider using tuple

// param pack only for size_t
template <size_t... Values>
struct UIntParamPack
{
    // PrependedBy<UIntParamPack<1, 2>>
    template <typename UIntParamPackToPrepend>
    struct PrependedBy
    {
    private:
        template <typename UIntParamPack2>
        struct _Base;
        template <size_t... Values2>
        struct _Base<UIntParamPack<Values2...>>
        {
            using Type = UIntParamPack<Values2..., Values...>;
        };

    public:
        using Type = typename _Base<UIntParamPackToPrepend>::Type;
    };

    // AppendedBy<UIntParamPack<3, 4>>
    template <typename UIntParamPackToAppend>
    struct AppendedBy
    {
    private:
        template <typename UIntParamPack2>
        struct _Base;
        template <size_t... Values2>
        struct _Base<UIntParamPack<Values2...>>
        {
            using Type = UIntParamPack<Values..., Values2...>;
        };
        
    public:
        using Type = typename _Base<UIntParamPackToAppend>::Type;
    };

    // At<0>
    template <size_t Index>
    struct At
    {
    private:
        template <size_t Index_, size_t Value_, size_t... Values_>
        struct _Base
        {
            static constexpr auto value = _Base<Index_ - 1, Values_...>::value;
        };
        
        template <size_t Value_, size_t... Values_>
        struct _Base<0, Value_, Values_...>
        {
            static constexpr auto value = Value_;
        };

    public:
        static constexpr auto value = _Base<Index, Values...>::value;
    };
};

// terminating UIntParamPack
template <>
struct UIntParamPack<>
{
    template <typename UIntParamPackToPrepend>
    struct PrependedBy
    {
        using Type = UIntParamPackToPrepend;
    };

    template <typename UIntParamPackToAppend>
    struct AppendedBy
    {
        using Type = UIntParamPackToAppend;
    };
};

struct Determinant
{
    template <size_t N>
    static float compute(const float (&m)[N][N])
    {
        return Base<N>::compute(m);
    }
    
    template <size_t M, size_t N>
    static float compute(const float (&m)[M][N])
    {
        static_assert( M == N, "requires square matrix" );
        return compute(m);
    }
    
    template <size_t Size>
    struct Base
    {
        // compute matrix determinant
        static float compute(const float (&m)[Size][Size])
        {
            using InitialElementIndices = typename InitialElementIndices<Size * Size>::Type;
            
            return Minor<Size, InitialElementIndices>::compute(m);
        }

        template <size_t M, size_t N>
        static float compute(const float (&m)[M][N])
        {
            ArrayValidator::perform(m);
            return compute(m);
        }
        
        struct ArrayValidator
        {
            template <size_t M, size_t N>
            static void perform(const float (&m)[M][N])
            {
                static_assert( M == Size && N == Size, "requires agreed dimensions" );
            }
        };

        // initial element indices, will be [0, 1, 2, 3, 4, 5, 6, 7, 8]
        // http://stackoverflow.com/a/10178791
        template <size_t N, typename Dummy = void>
        struct InitialElementIndices
        {
            using Type = typename InitialElementIndices<N - 1>::Type::template AppendedBy<UIntParamPack<N - 1>>::Type;
        };
        
        template <typename Dummy>
        struct InitialElementIndices<0, Dummy>
        {
            using Type = UIntParamPack<>;
        };

        // working with minors
        template <size_t MinorSize, typename ElementIndices>
        struct Minor
        {
            // compute minor determinant
            static float compute(const float (&m)[Size][Size])
            {
                return All<MinorSize>::compute(m);
            }

            template <size_t M, size_t N>
            static float compute(const float (&m)[M][N])
            {
                ArrayValidator::perform(m);
                return compute(m);
            }

            // building element indices for minor at (R, C)
            template <size_t R, size_t C, size_t N>
            struct MinorElementIndices
            {
                static constexpr auto CurrentIndex = MinorSize * MinorSize - N;
                static constexpr auto CurrentR = CurrentIndex / MinorSize;
                static constexpr auto CurrentC = CurrentIndex % MinorSize;
                
                using Rest = typename MinorElementIndices<R, C, N - 1>::Type;
                
                static constexpr auto ElementToSkip = (CurrentR != R && CurrentC != C);
                static constexpr auto CurrentElementIndex = ElementIndices::template At<CurrentIndex>::value;
                using Part = typename conditional<ElementToSkip, UIntParamPack<CurrentElementIndex>, UIntParamPack<>>::type;
                
                using Type = typename Rest::template AppendedBy<Part>::Type;
            };
            
            template <size_t R, size_t C>
            struct MinorElementIndices<R, C, 0>
            {
                using Type = UIntParamPack<>;
            };

            // loop all elements in the first row
            template <size_t Count, typename Dummy = void>
            struct All
            {
                static float compute(const float (&m)[Size][Size])
                {
                    constexpr auto CurrentIndex = Count - 1;
                    constexpr auto ElementIndex = ElementIndices::template At<CurrentIndex>::value;
                    
                    auto element = m[ElementIndex / Size][ElementIndex % Size];
                    auto minorDeterminant = Minor<MinorSize - 1, typename MinorElementIndices<0, CurrentIndex, MinorSize * MinorSize>::Type>::compute(m);
                    auto part = element * minorDeterminant;
                    
                    auto rest = All<Count - 1>::compute(m);
                    return (CurrentIndex % 2 == 0) ? (rest + part) : (rest - part);
                }

                template <size_t M, size_t N>
                static float compute(const float (&m)[M][N])
                {
                    ArrayValidator::perform(m);
                    return compute(m);
                }
            };
            
            template <typename Dummy>
            struct All<0, Dummy>
            {
                static float compute(const float (&m)[Size][Size])
                {
                    return 0;
                }
            };
        };
        
        template <typename ElementIndices>
        struct Minor<0, ElementIndices>
        {
            static float compute(const float (&m)[Size][Size])
            {
                return 1;
            }
        };
    };
};

struct Inverse
{
    template <size_t N>
    static void compute(const float (&src)[N][N], float (&dst)[N][N])
    {
        auto invDeterminant = 1 / Determinant::compute(src);
        Base<N>::template All<N * N>::compute(src, invDeterminant, dst);
    }
    
    template <size_t M, size_t N, size_t P, size_t Q>
    static float compute(const float (&src)[M][N], float (&dst)[P][Q])
    {
        static_assert( M == N, "requires square matrix" );
        static_assert( P == Q, "requires square matrix" );
        static_assert( M == P, "requires same dimensions" );
        return compute(src, dst);
    }

    template <size_t Size>
    struct Base
    {
        // loop all elements
        template <size_t Count, typename Dummy = void>
        struct All
        {
            static void compute(const float (&src)[Size][Size], float invDeterminant, float (&dst)[Size][Size])
            {
                constexpr auto CurrentIndex = Count - 1;
                constexpr auto R = CurrentIndex / Size;
                constexpr auto C = CurrentIndex % Size;
                
                using InitialElementIndices = typename Determinant::Base<Size>::template InitialElementIndices<Size * Size>::Type;
                using ElementIndices = typename Determinant::Base<Size>::template Minor<Size, InitialElementIndices>::template MinorElementIndices<R, C, Size * Size>::Type;
                //using ElementIndices = typename CofactorMatrixElementIndices<R, C, (Size - 1) * (Size - 1)>::Type;

                auto part = Determinant::Base<Size>::template Minor<Size - 1, ElementIndices>::compute(src);
                dst[C][R] = invDeterminant * (((R + C) % 2) == 0 ? 1 : -1) * part;                  // swapping R and C because we need the transpose of cofactor matrix
                All<Count - 1>::compute(src, invDeterminant, dst);
            }
            
            template <size_t M, size_t N, size_t P, size_t Q>
            static float compute(const float (&src)[M][N], float invDeterminant, float (&dst)[P][Q])
            {
                static_assert( M == N, "requires square matrix" );
                static_assert( P == Q, "requires square matrix" );
                static_assert( M == P, "requires same dimensions" );
                return compute(src, invDeterminant, dst);
            }
        };
        
        template <typename Dummy>
        struct All<0, Dummy>
        {
            static void compute(const float (&)[Size][Size], float, float (&)[Size][Size]) {}
        };

        // building element indices for cofactor matrix of m(R, C). reuse minor element indices but interchange first two columns
        template <size_t R, size_t C, size_t N>
        struct CofactorMatrixElementIndices
        {
            using InitialElementIndices = typename Determinant::Base<Size>::template InitialElementIndices<Size * Size>::Type;
            using MinorElementIndices = typename Determinant::Base<Size>::template Minor<Size, InitialElementIndices>::template MinorElementIndices<R, C, Size * Size>::Type;
            static constexpr auto CofactorMatrixSize = Size - 1;
            static constexpr auto CurrentIndex = CofactorMatrixSize * CofactorMatrixSize - N;

            static constexpr auto UnadjustedC = CurrentIndex % CofactorMatrixSize;

            static constexpr auto CurrentR = CurrentIndex / CofactorMatrixSize;
            static constexpr auto CurrentC = (UnadjustedC == 0 ? ((UnadjustedC + 1) % CofactorMatrixSize) : ((UnadjustedC == 1 ? (UnadjustedC - 1) : UnadjustedC)));
            
            static constexpr auto MappedIndex = CurrentR * CofactorMatrixSize + CurrentC;
            
            using Rest = typename CofactorMatrixElementIndices<R, C, N - 1>::Type;
            
            using Type = typename Rest::template AppendedBy<UIntParamPack<MappedIndex>>::Type;
        };
        
        template <size_t R, size_t C>
        struct CofactorMatrixElementIndices<R, C, 0>
        {
            using Type = UIntParamPack<>;
        };
    };
};

enum class TestEnum: size_t
{
    TestValue = UIntParamPack<1>::AppendedBy<UIntParamPack<2, 3>>::Type::AppendedBy<UIntParamPack<4>>::Type::At<2>::value,
};

int main(int argc, const char * argv[]) {
    
    static_assert( (size_t)TestEnum::TestValue == (size_t)3, "UIntParamPack should work properly" );

    float m[2][2] = {
        {1, 2},
        {3, 4},
    };
//    float m2[3][3] = {
//        {1, 3, 5},
//        {7, 9, 1},
//        {3, 5, 7},
//    };
    float m2[3][3] = {
        {1, 3, 5},
        {7, 9, 1},
        {3, 5, 7},
    };
    float m3[4][4] = {
        {2, 0, 0, 0},
        {0, 2, 0, 0},
        {0, 0, 2, 0},
        {1, 2, 3, 1},
    };
    
    cout << Determinant::compute(m) << endl;
    cout << Determinant::compute(m2) << endl;
    cout << Determinant::compute(m3) << endl;
    
    decltype(m) m4;
    decltype(m2) m5;
    decltype(m3) m6;

    Inverse::compute(m, m4);
    Inverse::compute(m2, m5);
    Inverse::compute(m3, m6);

    return 0;
}
