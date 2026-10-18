// Force-included when the native tests compile the engine-independent sources. Unreal and the
// Windows headers it pulls in define these as macros; if any of our identifiers collide with
// one, the native build fails here instead of inside Unreal Build Tool.
#pragma once
#define PI (3.1415926535897932f)
#define SMALL_NUMBER (1.e-8f)
#define KINDA_SMALL_NUMBER (1.e-4f)
#define BIG_NUMBER (3.4e+38f)
#define DELTA (0.00001f)
#define INDEX_NONE (-1)
#define IN
#define OUT
#define TEXT(x) L##x
#define check(expr) do { (void)(expr); } while (0)
#define verify(expr) do { (void)(expr); } while (0)
#define ensure(expr) (!!(expr))
#define CONSTEXPR constexpr
#define FORCEINLINE inline
#define UE_PI (3.14159265358979323846)
#define HALF_PI (1.57079632679f)
#define TWO_PI (6.28318530718f)
#define INV_PI (0.31830988618f)
#define EULERS_NUMBER (2.71828182845904523536f)
#define near
#define far
#define small char
#define hyper __int64
#define ERROR 0
#define DELETE 0x00010000L
#define IGNORE 0
#define INFINITE 0xFFFFFFFF
#define OPTIONAL
#define CONST const
#define VOID void
#define ABSOLUTE 1
#define RELATIVE 2
#define TRANSPARENT 1
#define OPAQUE 2
#define DIFFERENCE 11
#define PASCAL
#define CALLBACK
#define interface struct
#define TRUE 1
#define FALSE 0
#define LIKELY(x) (x)
#define UNLIKELY(x) (x)
#define MAX_int32 2147483647
#define MAX_flt (3.402823466e+38F)
