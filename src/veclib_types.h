/******************************************************************************/
/*                                                                            */
/* Licensed Materials - Property of IBM                                       */
/*                                                                            */
/* IBM Power Vector Intrinisic Functions version 1.0.4                        */
/*                                                                            */
/* Copyright IBM Corp. 2014,2016                                              */
/* US Government Users Restricted Rights - Use, duplication or                */
/* disclosure restricted by GSA ADP Schedule Contract with IBM Corp.          */
/*                                                                            */
/* See the licence in the license subdirectory.                               */
/*                                                                            */
/* More information on this software is available on the IBM DeveloperWorks   */
/* website at                                                                 */
/*  https://www.ibm.com/developerworks/community/groups/community/powerveclib */
/*                                                                            */
/******************************************************************************/

#ifndef _H_VECLIB_TYPES
#define _H_VECLIB_TYPES

/******************************************************************************/

/** Control Macros **/

/* Control alignment */
#define VECLIB_ALIGNED8 __attribute__ ((__aligned__ (8)))
#define VECLIB_ALIGNED16 __attribute__ ((__aligned__ (16)))
/* Note: PowerPC only needs 16 byte alignment and stacks are only aligned 16 */
/* On some compilers: __declspec (align (N)) */
/* In C++ 11:         alignas (N) */

/* Control inlining */
#ifdef NOT_ALWAYS_INLINE
  #define VECLIB_INLINE
#else
  #define VECLIB_INLINE static inline __attribute__ ((__always_inline__))
#endif
#define VECLIB_NOINLINE static __attribute__ ((__noinline__))

/* Control compiler */
#if defined(__xlc__) || defined(__xlC__)
/* AIX xlc defines __IBMC__ or __xlc__; OpenPOWER defines __ibmxl__ */
/* AIX xlC defines __IBMCPP__ or __xlC__; OpenPOWER defines __ibmxl__ */
  #ifndef __ibmxl__
    #define __ibmxl__
  #endif
#endif

#ifdef __GNUC__
  /* GCC version */
  #define __GCC_VERSION__ (__GNUC__ * 100 + __GNUC_MINOR__ * 10 + __GNUC_PATCHLEVEL__)
#endif

/* Control double precision */
#ifdef _ARCH_PWR7
  #define VECLIB_VSX
#endif

/* Control BIFs */

/* Gather bits by bytes by doubleword */
#if defined (__ibmxl__)
  #if defined (__LITTLE_ENDIAN__)
    #define VEC_GB(v) ( vec_gb ((vector unsigned char) v) )
  #else
    #define VEC_GB(v) ( (vector unsigned char) vec_gbb ((vector unsigned long long) v) )
  #endif
#elif defined (__clang__)
  #define VEC_GB(v) vec_gb (v)
#elif defined (__GNUC__)
  #if __GCC_VERSION__ >= 600
    #define VEC_GB(v) vec_gb ((vector unsigned char) v)
  #else
    #define VEC_GB(v) vec_vgbbd ((vector unsigned char) v)
  #endif
#else
  #error Compiler not supported yet.
#endif

/* Count leading zeros */
#if defined (__ibmxl__)
  #define VEC_CNTLZ vec_cntlz
#elif defined (__clang__)
  #define VEC_CNTLZ vec_cntlz
#elif defined (__GNUC__)
  #if __GCC_VERSION__ >= 501
    #define VEC_CNTLZ vec_cntlz
  #else
    #define VEC_CNTLZ vec_vclz
  #endif
#else
  #error Compiler not supported yet.
#endif


/** Types **/

typedef
  VECLIB_ALIGNED8
  unsigned long long
__m64;

typedef
  VECLIB_ALIGNED16
  vector float
__m128;

#ifdef VECLIB_VSX
  typedef
    VECLIB_ALIGNED16
    vector double
  __m128d;
#endif

typedef
  VECLIB_ALIGNED16
  vector unsigned char
__m128i;

typedef
  VECLIB_ALIGNED16
  struct {
    __m128i m128i_0;
    __m128i m128i_1;
  } __m256i;

typedef
  VECLIB_ALIGNED16
  struct {
    __m128 m128_0;
    __m128 m128_1;
  } __m256;

#ifdef VECLIB_VSX
  typedef
    VECLIB_ALIGNED16
    struct {
      __m128d m128d_0;
      __m128d m128d_1;
    } __m256d;
#endif


/** Parameter Types **/

typedef const long intlit;   /* general int literal */
typedef const long intlit1;  /* 1 bit int literal */
typedef const long intlit2;  /* 2 bit int literal */
typedef const long intlit3;  /* 3 bit int literal */
typedef const long intlit4;  /* 4 bit int literal */
typedef const long intlit5;  /* 5 bit int literal */
typedef const long intlit8;  /* 8 bit int literal */

typedef const long uintlit;   /* general unsigned int literal */
typedef const long uintlit1;  /* 1 bit unsigned int literal */
typedef const long uintlit2;  /* 2 bit unsigned int literal */
typedef const long uintlit3;  /* 3 bit unsigned int literal */
typedef const long uintlit4;  /* 4 bit unsigned int literal */
typedef const long uintlit5;  /* 5 bit unsigned int literal */
typedef const long uintlit8;  /* 8 bit unsigned int literal */


/** Internal Types and Constants **/

typedef
  VECLIB_ALIGNED8
  union {
    __m64                     as_m64;
    char                      as_char              [8];
    signed char               as_signed_char       [8];
    short                     as_short             [4];
    int                       as_int               [2];
    long long                 as_long_long;
    float                     as_float             [2];
    double                    as_double;
  } __m64_union;

typedef
  VECLIB_ALIGNED8
  union {
    __m64                     as_m64;
    signed char               as_char              [8];
    unsigned char             as_unsigned_char     [8];
    signed char               as_signed_char       [8];
    unsigned short            as_unsigned_short    [4];
    short                     as_short             [4];
    int                       as_int               [2];
    long long                 as_long_long;
    float                     as_float             [2];
    double                    as_double;
  } __m64_all_union;

typedef
  VECLIB_ALIGNED16
  union {
    __m128i                   as_m128i;
    __m64                     as_m64               [2];
    vector signed   char      as_vector_signed_char;
    vector unsigned char      as_vector_unsigned_char;
    vector bool     char      as_vector_bool_char;
    vector signed   short     as_vector_signed_short;
    vector unsigned short     as_vector_unsigned_short;
    vector bool     short     as_vector_bool_short;
    vector signed   int       as_vector_signed_int;
    vector unsigned int       as_vector_unsigned_int;
    vector bool     int       as_vector_bool_int;
    vector signed   long long as_vector_signed_long_long;
    vector unsigned long long as_vector_unsigned_long_long;
    vector bool     long long as_vector_bool_long_long;
    char                      as_char              [16];
    short                     as_short             [8];
    int                       as_int               [4];
    unsigned int              as_unsigned_int      [4];
    long long                 as_long_long         [2];
  } __m128i_union;

typedef
  VECLIB_ALIGNED16
  union {
    __m128                    as_m128;
    __m64                     as_m64               [2];
    vector float              as_vector_float;
    float                     as_float             [4];
    int                       as_hex               [4];
  } __m128_union;

typedef
  VECLIB_ALIGNED16
  union {
    __m128i                   as_m128i;
    __m64                     as_m64               [2];
    vector signed   char      as_vector_signed_char;
    vector unsigned char      as_vector_unsigned_char;
    vector bool     char      as_vector_bool_char;
    vector signed   short     as_vector_signed_short;
    vector unsigned short     as_vector_unsigned_short;
    vector bool     short     as_vector_bool_short;
    vector signed   int       as_vector_signed_int;
    vector unsigned int       as_vector_unsigned_int;
    vector bool     int       as_vector_bool_int;
    vector signed   long long as_vector_unsigned_long_long;
    vector unsigned long long as_vector_signed_long_long;
    vector bool     long long as_vector_bool_long_long;
    char                      as_char              [16];
    short                     as_short             [8];
    int                       as_int               [4];
    long long                 as_long_long         [2];
    __m128                    as_m128;
    vector float              as_vector_float;
    float                     as_float             [4];
    #ifdef VECLIB_VSX
      __m128d                 as_m128d;
      vector double           as_vector_double;
      double                  as_double            [4];
    #endif
  } __m128_all_union;

#ifdef VECLIB_VSX
  typedef
    VECLIB_ALIGNED16
    union {
      __m128d                 as_m128d;
      __m64                   as_m64               [2];
      vector double           as_vector_double;
      double                  as_double            [2];
      int                     as_int               [4];
      long long               as_long_hex          [2];
    } __m128d_union;
#endif

typedef
  VECLIB_ALIGNED16
  union {
    __m256i                   as_m256i;
    __m128i                   as_m128i                      [4];
    __m64                     as_m64                        [4];
    vector signed   char      as_vector_signed_char         [2];
    vector unsigned char      as_vector_unsigned_char       [2];
    vector bool     char      as_vector_bool_char           [2];
    vector signed   short     as_vector_signed_short        [2];
    vector unsigned short     as_vector_unsigned_short      [2];
    vector bool     short     as_vector_bool_short          [2];
    vector signed   int       as_vector_signed_int          [2];
    vector unsigned int       as_vector_unsigned_int        [2];
    vector bool     int       as_vector_bool_int            [2];
    vector signed   long long as_vector_signed_long_long    [2];
    vector unsigned long long as_vector_unsigned_long_long  [2];
    vector bool     long long as_vector_bool_long_long      [2];
    char                      as_char                       [32];
    short                     as_short                      [16];
    int                       as_int                        [8];
    unsigned int              as_unsigned_int               [8];
    long long                 as_long_long                  [4];
  } __m256i_union;

typedef
  VECLIB_ALIGNED16
  union {
    __m256                    as_m256;
    __m128                    as_m128              [2];
    __m64                     as_m64               [4];
    vector float              as_vector_float      [2];
    float                     as_float             [8];
    int                       as_hex               [8];
  } __m256_union;

typedef
  VECLIB_ALIGNED16
  union {
    __m256i                   as_m256i;
    __m256                    as_m256;
    __m128d                   as_m128d                      [2];
    __m128                    as_m128                       [2];
    __m128i                   as_m128i                      [2];
    __m64                     as_m64                        [4];
    vector signed   char      as_vector_signed_char         [2];
    vector unsigned char      as_vector_unsigned_char       [2];
    vector bool     char      as_vector_bool_char           [2];
    vector signed   short     as_vector_signed_short        [2];
    vector unsigned short     as_vector_unsigned_short      [2];
    vector bool     short     as_vector_bool_short          [2];
    vector signed   int       as_vector_signed_int          [2];
    vector unsigned int       as_vector_unsigned_int        [2];
    vector bool     int       as_vector_bool_int            [2];
    vector signed   long long as_vector_unsigned_long_long  [2];
    vector unsigned long long as_vector_signed_long_long    [2];
    vector bool     long long as_vector_bool_long_long      [2];
    vector float              as_vector_float               [2];
    char                      as_char                       [32];
    short                     as_short                      [16];
    int                       as_int                        [8];
    float                     as_float                      [8];
    long long                 as_long_long                  [4];
    #ifdef VECLIB_VSX
      double                  as_double                     [4];
      vector double           as_vector_double              [2];
      __m256d                 as_m256d;
    #endif
  } __m256_all_union;

#ifdef VECLIB_VSX
  typedef
    VECLIB_ALIGNED16
    union {
      __m256d                 as_m256d;
      __m128d                 as_m128d             [2];
      __m64                   as_m64               [4];
      vector double           as_vector_double     [2];
      double                  as_double            [4];
      long long               as_long_hex          [4];
    } __m256d_union;
#endif

/** Element Position Macros **/

#define __vsr_lowest_left_half_char  7
#define __vsr_lowest_left_half_short 3
#define __vsr_lowest_left_half_int   1
#define __vsr_left_half_long_long    0

#ifdef __LITTLE_ENDIAN__
  #define __vsx_lowest_left_half_char_in_memory  8
  #define __vsr_lowest_left_half_short_in_memory 4
  #define __vsr_lowest_left_half_int_in_memory   2
  #define __vsr_left_half_long_long_in_memory    1
  #define __vsr_right_half_long_long_in_memory   0
#elif __BIG_ENDIAN__
  #define __vsx_lowest_left_half_char_in_memory  7
  #define __vsr_lowest_left_half_short_in_memory 3
  #define __vsr_lowest_left_half_int_in_memory   1
  #define __vsr_left_half_long_long_in_memory    0
  #define __vsr_right_half_long_long_in_memory   1
#endif

/******************************************************* Permute ******************************************************/

#define _MM_SHUFFLE(a, b, c, d) ((intlit8) (a*64+b*16+c*4+d))

#endif
