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

#ifndef _H_VEC128INT
#define _H_VEC128INT

#include <altivec.h>
#include "veclib_types.h"

/******************************************************** Load ********************************************************/

/* Load 128-bits of integer data, aligned */
VECLIB_INLINE __m128i vec_load1q (__m128i const* address)
{
  return (__m128i) vec_ld (0, (vector unsigned char*) address);
}

/* Load 128-bits of integer data, unaligned */
VECLIB_INLINE __m128i vec_loadu1q (__m128i const* address)
{
  __m128i src1, src2;
  src1 = vec_ld(0, (vector unsigned char*) address);
  src2 = vec_ld(16, (vector unsigned char*) address);

  #if (defined(__ibmxl__) && (defined(__LITTLE_ENDIAN__) || defined(__BIG_ENDIAN__)))
  return (__m128i) vec_perm (src1, src2, vec_lvsl(0,(unsigned char*)address));
  #elif ((defined __GNUC__) && (__GCC_VERSION__ >= 492))
    return vec_xl (0, (unsigned char *) address);
  #elif (defined __GNUC__) && (__GCC_VERSION__ < 492)
    static const vector unsigned char permMask = { 0x0F, 0x0E, 0x0D, 0x0C, 0x0B, 0x0A, 0x09, 0x08, 0x07, 0x06, 0x05, 0x04, 0x03, 0x02, 0x01, 0x00 };
    return (__m128i) vec_perm (src1, src2, vec_perm( vec_lvsl(0, ( unsigned char*)address), src1, permMask));
  #endif
}

/* Load 128-bits of integer data, unaligned - deprecated - use previous function */
VECLIB_INLINE __m128i vec_load1qu (__m128i const* address)
{
  return vec_loadu1q (address);
}

/* Load 64-bits of integer data to lower part and zero upper part */
VECLIB_INLINE __m128i vec_loadlower1sd (__m128i const* from) {
  #if (defined(__ibmxl__) && defined(__LITTLE_ENDIAN__))
    vector bool long long mask = (vector bool long long) { 0xFFFFFFFFFFFFFFFF, 0x0000000000000000};
    return (__m128i)vec_and( vec_loadu1q(from) , (vector unsigned char)mask);
  #elif (defined(__ibmxl__) && defined(__BIG_ENDIAN__))
    vector bool long long mask = (vector bool long long) { 0x0000000000000000, 0xFFFFFFFFFFFFFFFF};
    return (__m128i)vec_and(*from, (vector unsigned char)mask);
  #elif defined(__GNUC__)
    vector bool long long mask = (vector bool long long) { 0xFFFFFFFFFFFFFFFF, 0x0000000000000000};
    return (__m128i)vec_and(*from, (vector unsigned char)mask);
  #endif
}

/******************************************************** Set *********************************************************/

/* Set 128 integer bits to zero */
VECLIB_INLINE __m128i vec_zero1q (void)
{
  return (__m128i) (vector signed int) vec_splats (0);
}

/* Splat 8-bit char to 16 8-bit chars */
VECLIB_INLINE __m128i vec_splat16sb (char scalar)
{ return (__m128i) vec_splats ((signed char) scalar); }

/* Splat 16-bit short to 8 16-bit shorts */
VECLIB_INLINE __m128i vec_splat8sh (short scalar)
{ return (__m128i) vec_splats (scalar); }

/* Splat 32-bit int to 4 32-bit ints */
VECLIB_INLINE __m128i vec_splat4sw (int scalar)
{ return (__m128i) vec_splats (scalar); }

/* Splat 64-bit long long to 2 64-bit long longs */
VECLIB_INLINE __m128i vec_splat2sd (long long scalar)
{ return (__m128i) vec_splats (scalar); }

/* Set 16 8-bit chars */
VECLIB_INLINE __m128i vec_set16sb (char c15, char c14, char c13, char c12, char c11, char c10, char c9, char c8, char c7, char c6, char c5, char c4, char c3, char c2, char c1, char c0)
{
  __m128i_union t;
  #if __LITTLE_ENDIAN__
    t.as_char[0] = c0;
    t.as_char[1] = c1;
    t.as_char[2] = c2;
    t.as_char[3] = c3;
    t.as_char[4] = c4;
    t.as_char[5] = c5;
    t.as_char[6] = c6;
    t.as_char[7] = c7;
    t.as_char[8] = c8;
    t.as_char[9] = c9;
    t.as_char[10] = c10;
    t.as_char[11] = c11;
    t.as_char[12] = c12;
    t.as_char[13] = c13;
    t.as_char[14] = c14;
    t.as_char[15] = c15;
  #elif __BIG_ENDIAN__
    t.as_char[0] = c15;
    t.as_char[1] = c14;
    t.as_char[2] = c13;
    t.as_char[3] = c12;
    t.as_char[4] = c11;
    t.as_char[5] = c10;
    t.as_char[6] = c9;
    t.as_char[7] = c8;
    t.as_char[8] = c7;
    t.as_char[9] = c6;
    t.as_char[10] = c5;
    t.as_char[11] = c4;
    t.as_char[12] = c3;
    t.as_char[13] = c2;
    t.as_char[14] = c1;
    t.as_char[15] = c0;
  #endif
  return (__m128i) t.as_vector_unsigned_char;
}

/* Set 8 16-bit shorts */
VECLIB_INLINE __m128i vec_set8sh (short s7, short s6, short s5, short s4, short s3, short s2, short s1, short s0)
{
  __m128i_union t;
  #if __LITTLE_ENDIAN__
    t.as_short[0] = s0;
    t.as_short[1] = s1;
    t.as_short[2] = s2;
    t.as_short[3] = s3;
    t.as_short[4] = s4;
    t.as_short[5] = s5;
    t.as_short[6] = s6;
    t.as_short[7] = s7;
  #elif __BIG_ENDIAN__
    t.as_short[0] = s7;
    t.as_short[1] = s6;
    t.as_short[2] = s5;
    t.as_short[3] = s4;
    t.as_short[4] = s3;
    t.as_short[5] = s2;
    t.as_short[6] = s1;
    t.as_short[7] = s0;
  #endif
  return (__m128i) t.as_vector_unsigned_short;
}

/* Set 4 32-bit ints */
VECLIB_INLINE __m128i vec_set4sw (int i3, int i2, int i1, int i0)
{
  __m128i_union t;
  #ifdef __LITTLE_ENDIAN__
    t.as_int[0] = i0;
    t.as_int[1] = i1;
    t.as_int[2] = i2;
    t.as_int[3] = i3;
  #elif __BIG_ENDIAN__
    t.as_int[0] = i3;
    t.as_int[1] = i2;
    t.as_int[2] = i1;
    t.as_int[3] = i0;
  #endif
  return (__m128i) t.as_vector_signed_int;
}

/* Set 2 64-bit long longs */
VECLIB_INLINE __m128i vec_set2sd (__m64 LL1, __m64 LL0)
{
  __m64_union LL0_union;
  __m64_union LL1_union;
  LL0_union.as_m64 = LL0;
  LL1_union.as_m64 = LL1;
  __m128i_union t;
  #ifdef __LITTLE_ENDIAN__
    t.as_long_long[0] = LL0_union.as_long_long;
    t.as_long_long[1] = LL1_union.as_long_long;
  #elif __BIG_ENDIAN__
    t.as_long_long[0] = LL1_union.as_long_long;
    t.as_long_long[1] = LL0_union.as_long_long;
  #endif
  /* vector signed long long is being phased in for gcc */
  return (__m128i) t.as_vector_unsigned_char;
}

/* Set 16 8-bit chars reversed */
VECLIB_INLINE __m128i vec_setreverse16sb (char c15, char c14, char c13, char c12, char c11, char c10, char c9, char c8, char c7, char c6, char c5, char c4, char c3, char c2, char c1, char c0)
{
  return vec_set16sb (c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14, c15);
}

/* Set 8 16-bit shorts reversed */
VECLIB_INLINE __m128i vec_setreverse8sh (short s7, short s6, short s5, short s4, short s3, short s2, short s1, short s0)
{
  return vec_set8sh (s0, s1, s2, s3, s4, s5, s6, s7);
}

/* Set 4 32-bit ints reversed */
VECLIB_INLINE __m128i vec_setreverse4sw (int i3, int i2, int i1, int i0)
{
  return vec_set4sw (i0, i1, i2, i3);
}

/* Set 2 64-bit long longs reversed */
VECLIB_INLINE __m128i vec_setreverse2sd (__m64 v1, __m64 v0)
{
  return vec_set2sd (v0, v1);
}

/* Set lower 64-bits of integer data and zero upper part */
VECLIB_INLINE __m128i vec_Zerouppersd (__m128i v) {
  #ifdef __LITTLE_ENDIAN__
    #if ( defined __ibmxl__ ) || (defined __GNUC__) && (__GCC_VERSION__ >= 492)
      vector signed int zeroes = {0x00000000, 0x00000000, 0x00000000, 0x00000000 };
      return (__m128i)vec_mergeh((vector unsigned long long)v, (vector unsigned long long)zeroes);
    #elif (defined __GNUC__) && (__GCC_VERSION__ < 492)
      vector signed int mask = {0xffffffff, 0xffffffff, 0x00000000, 0x00000000 };
      return (__m128i)vec_and((vector unsigned int)v, (vector unsigned int)mask);
    #endif
  #elif __BIG_ENDIAN__
    vector signed int mask = {0x00000000, 0x00000000, 0xffffffff, 0xffffffff };
    return (__m128i)vec_and( (vector unsigned long long)v, (vector unsigned long long)mask );
  #endif
}

/******************************************************** Store *******************************************************/

/* Store 128-bits integer, aligned */
VECLIB_INLINE void vec_store1q (__m128i* address, __m128i v)
{ vec_st (v, 0, address); }

/* Store 128-bits integer, unaligned */
VECLIB_INLINE void vec_storeu1q (__m128i* to, __m128i from)
{
  #ifdef __LITTLE_ENDIAN__
    #ifdef __ibmxl__
      vec_xst (from, 0, (unsigned char*) to);
    #elif (defined __GNUC__) && (__GCC_VERSION__ >= 492)
      vec_xst (from, 0, (unsigned char*) to);
    #elif (defined __GNUC__) && (__GCC_VERSION__ < 492)
      __m128_all_union from_union; from_union.as_m128i = from;

      /* Prepare for later generate select control mask vector */
      vector signed char all_one = vec_splat_s8( -1 );
      vector signed char all_zero = vec_splat_s8( 0 );
      /* Generate permute vector for the upper part of each m128d component */
      vector unsigned char permute_vector = vec_lvsr (0, (unsigned char *) to);
      /* Generate selector vector for the lower part of each m128d component */
      vector unsigned char select_vector = vec_perm ((vector unsigned char) all_zero, (vector unsigned char) all_one, permute_vector);

      /* Load from */
      /* Perform a 16-byte load of the original data starting at BoundAlign (to + 0) and BoundAlign (to + 16)*/
      vector unsigned char low = vec_ld (0, (unsigned char *) to);
      vector unsigned char high = vec_ld (16, (unsigned char *) to);
      /* Perform permute, the result will look like:
         original data ... from_union.as_vector_unsigned_char ... original data */
      vector unsigned char temp_low = vec_perm (from_union.as_vector_unsigned_char, from_union.as_vector_unsigned_char, permute_vector);
      low = vec_sel (low, temp_low, select_vector);
      high = vec_perm (from_union.as_vector_unsigned_char, high, permute_vector);
      /* Store the aligned result for from_union.as_m128 */
      vec_st (low, 0, (unsigned char *) to);
      vec_st (high, 16, (unsigned char *) to);
    #endif
  #elif __BIG_ENDIAN__
    /* Prepare for later generate control mask vector */
    vector signed char all_one = vec_splat_s8( -1 );
    vector signed char all_zero = vec_splat_s8( 0 );
    /* Generate permute vector for the upper part of from component */
    vector unsigned char permute_vector = vec_lvsr (0, (unsigned char *) to);
    /* Generate selector vector for the lower part of from component */
    vector unsigned char select_vector = vec_perm ((vector unsigned char) all_zero, (vector unsigned char) all_one, permute_vector);

    /* Perform a 16-byte load of the original data starting at BoundAlign (to + 0) and BoundAlign (to + 16)*/
    vector unsigned char low = vec_ld (0, (unsigned char *) to);
    vector unsigned char high = vec_ld (16, (unsigned char *) to);
    /* Perform permute, the result will look like: original data ... from ... original data */
    vector unsigned char temp_low = vec_perm (low, from, permute_vector);
    low = vec_sel (low, temp_low, select_vector);
    high = vec_perm (from, high, permute_vector);
    /* Store the aligned result for from */
    vec_st (low, 0, (unsigned char *) to);
    vec_st (high, 16, (unsigned char *) to);
  #endif
}

/* Store 128-bits integer, unaligned - deprecated - use previous function */
VECLIB_INLINE void vec_store1qu (__m128i* to, __m128i from)
{
  vec_storeu1q (to, from);
}

/* Store lower 64-bit long long */
VECLIB_INLINE void vec_storelower1sdof2sd (__m128i* to, __m128i from)
{
  __m128i_union from_union; from_union.as_m128i = from;
  __m128i_union to_union; to_union.as_m128i = *to;
  __m128i_union ret_union;
  #ifdef __LITTLE_ENDIAN__
    ret_union.as_m64[0] = from_union.as_m64[0];
    ret_union.as_m64[1] = to_union.as_m64[1];
  #elif __BIG_ENDIAN__
    ret_union.as_m64[1] = from_union.as_m64[1];
    ret_union.as_m64[0] = to_union.as_m64[0];
  #endif
  *to = ret_union.as_m128i;
}

/******************************************************* Insert *******************************************************/

/* Insert 32-bit int */
VECLIB_INLINE __m128i vec_insert4sw (__m128i into, int from, const intlit2 element_from_right)
{
  static const vector unsigned char permute_selector[4] = {
    /* To insert lowest element into specified element with other elements unchanged */
    #ifdef __LITTLE_ENDIAN__
      { 0x10,0x11,0x12,0x13, 0x04,0x05,0x06,0x07, 0x08,0x09,0x0A,0x0B, 0x0C,0x0D,0x0E,0x0F }, /* element 0 */
      { 0x00,0x01,0x02,0x03, 0x14,0x15,0x16,0x17, 0x08,0x09,0x0A,0x0B, 0x0C,0x0D,0x0E,0x0F }, /* element 1 */
      { 0x00,0x01,0x02,0x03, 0x04,0x05,0x06,0x07, 0x18,0x19,0x1A,0x1B, 0x0C,0x0D,0x0E,0x0F }, /* element 2 */
      { 0x00,0x01,0x1E,0x1F, 0x04,0x05,0x06,0x07, 0x08,0x09,0x0A,0x0B, 0x1C,0x1D,0x1E,0x1F }  /* element 3 */
    #elif __BIG_ENDIAN__
      { 0x00,0x01,0x02,0x03, 0x04,0x05,0x06,0x07, 0x08,0x09,0x0A,0x0B, 0x1C,0x1D,0x1E,0x1F }, /* element 0 */
      { 0x00,0x01,0x02,0x03, 0x04,0x05,0x06,0x07, 0x1C,0x1D,0x1E,0x1F, 0x0C,0x0D,0x0E,0x0F }, /* element 1 */
      { 0x00,0x01,0x02,0x03, 0x1C,0x1D,0x1E,0x1F, 0x08,0x09,0x0A,0x0B, 0x0C,0x0D,0x0E,0x0F }, /* element 2 */
      { 0x1C,0x1D,0x1E,0x1F, 0x04,0x05,0x06,0x07, 0x08,0x09,0x0A,0x0B, 0x0C,0x0D,0x0E,0x0F }  /* element 3 */
    #endif
  };
  return (__m128i) vec_perm ((vector signed int) into, vec_splats (from), permute_selector[element_from_right]);
}

/* Insert 64-bit long long */
VECLIB_INLINE __m128i vec_insert2sd (__m128i into, long long from, const intlit1 element_from_right)
{
  static const vector unsigned char permute_selector[2] = {
    /* To insert lowest element into specified element with other elements unchanged */
    #ifdef __LITTLE_ENDIAN__
      { 0x10,0x11,0x12,0x13, 0x14,0x15,0x16,0x17, 0x08,0x09,0x0A,0x0B, 0x0C,0x0D,0x0E,0x0F }, /* element 0 */
      { 0x00,0x01,0x02,0x03, 0x04,0x05,0x06,0x07, 0x18,0x19,0x1A,0x1B, 0x1C,0x1D,0x1E,0x1F }  /* element 1 */
    #elif __BIG_ENDIAN__
      { 0x00,0x01,0x02,0x03, 0x04,0x05,0x06,0x07, 0x18,0x19,0x1A,0x1B, 0x1C,0x1D,0x1E,0x1F }, /* element 0 */
      { 0x18,0x19,0x1A,0x1B, 0x1C,0x1D,0x1E,0x1F, 0x08,0x09,0x0A,0x0B, 0x0C,0x0D,0x0E,0x0F }  /* element 1 */
    #endif
  };
  __m128i_union from_union;
  from_union.as_long_long[0] = from;
  from_union.as_long_long[1] = from;
  return (__m128i) vec_perm (into, from_union.as_m128i, permute_selector[element_from_right]);
}

/* Insert 32-bit int, zeroing upper */
VECLIB_INLINE __m128i vec_convert1swto1uq (int from)
{ return vec_set4sw (0, 0, 0, from); }

/* Insert 16-bit short into one of 8 16-bit shorts */
VECLIB_INLINE __m128i vec_insert8sh (__m128i v, int scalar, intlit3 element_from_right)
{
  static const vector unsigned char permute_selector[8] = {
  /* To insert halfword from lowest halfword of the left half into specified element with other elements unchanged */
    #ifdef __LITTLE_ENDIAN__
      { 0x18,0x19, 0x02,0x03, 0x04,0x05, 0x06,0x07, 0x08,0x09, 0x0A,0x0B, 0x0C,0x0D, 0x0E,0x0F },  /* element 0 */
      { 0x00,0x01, 0x18,0x19, 0x04,0x05, 0x06,0x07, 0x08,0x09, 0x0A,0x0B, 0x0C,0x0D, 0x0E,0x0F },  /* element 1 */
      { 0x00,0x01, 0x02,0x03, 0x18,0x19, 0x06,0x07, 0x08,0x09, 0x0A,0x0B, 0x0C,0x0D, 0x0E,0x0F },  /* element 2 */
      { 0x00,0x01, 0x02,0x03, 0x04,0x05, 0x18,0x19, 0x08,0x09, 0x0A,0x0B, 0x0C,0x0D, 0x0E,0x0F },  /* element 3 */
      { 0x00,0x01, 0x02,0x03, 0x04,0x05, 0x06,0x07, 0x18,0x19, 0x0A,0x0B, 0x0C,0x0D, 0x0E,0x0F },  /* element 4 */
      { 0x00,0x01, 0x02,0x03, 0x04,0x05, 0x06,0x07, 0x08,0x09, 0x18,0x19, 0x0C,0x0D, 0x0E,0x0F },  /* element 5 */
      { 0x00,0x01, 0x02,0x03, 0x04,0x05, 0x06,0x07, 0x08,0x09, 0x0A,0x0B, 0x18,0x19, 0x0E,0x0F },  /* element 6 */
      { 0x00,0x01, 0x02,0x03, 0x04,0x05, 0x06,0x07, 0x08,0x09, 0x0A,0x0B, 0x0C,0x0D, 0x18,0x19 }   /* element 7 */
    #elif __BIG_ENDIAN__
      { 0x00,0x01, 0x02,0x03, 0x04,0x05, 0x06,0x07, 0x08,0x09, 0x0A,0x0B, 0x0C,0x0D, 0x16,0x17 },  /* element 0 */
      { 0x00,0x01, 0x02,0x03, 0x04,0x05, 0x06,0x07, 0x08,0x09, 0x0A,0x0B, 0x16,0x17, 0x0E,0x0F },  /* element 1 */
      { 0x00,0x01, 0x02,0x03, 0x04,0x05, 0x06,0x07, 0x08,0x09, 0x16,0x17, 0x0C,0x0D, 0x0E,0x0F },  /* element 2 */
      { 0x00,0x01, 0x02,0x03, 0x04,0x05, 0x06,0x07, 0x16,0x17, 0x0A,0x0B, 0x0C,0x0D, 0x0E,0x0F },  /* element 3 */
      { 0x00,0x01, 0x02,0x03, 0x04,0x05, 0x16,0x17, 0x08,0x09, 0x0A,0x0B, 0x0C,0x0D, 0x0E,0x0F },  /* element 4 */
      { 0x00,0x01, 0x02,0x03, 0x16,0x17, 0x06,0x07, 0x08,0x09, 0x0A,0x0B, 0x0C,0x0D, 0x0E,0x0F },  /* element 5 */
      { 0x00,0x01, 0x16,0x17, 0x04,0x05, 0x06,0x07, 0x08,0x09, 0x0A,0x0B, 0x0C,0x0D, 0x0E,0x0F },  /* element 6 */
      { 0x16,0x17, 0x02,0x03, 0x04,0x05, 0x06,0x07, 0x08,0x09, 0x0A,0x0B, 0x0C,0x0D, 0x0E,0x0F }   /* element 7 */
    #endif
  };
  __m128i_union t;
  t.as_long_long[__vsr_left_half_long_long_in_memory] = scalar;
  return (__m128i) vec_perm (v, t.as_m128i, permute_selector[element_from_right]);
}

/* Insert 8-bit unsigned char into one of 16 bytes */
VECLIB_INLINE __m128i vec_insert16ub (__m128i v, int scalar, intlit4 element_from_right)
{
  static const vector unsigned char permute_selector[16] = {
  /* To insert halfword from lowest halfword of the left half into specified element with other elements unchanged */
    #ifdef __LITTLE_ENDIAN__
      { 0x18,0x01, 0x02,0x03, 0x04,0x05, 0x06,0x07, 0x08,0x09, 0x0A,0x0B, 0x0C,0x0D, 0x0E,0x0F },  /* element 0 */
      { 0x00,0x18, 0x02,0x03, 0x04,0x05, 0x06,0x07, 0x08,0x09, 0x0A,0x0B, 0x0C,0x0D, 0x0E,0x0F },  /* element 1 */
      { 0x00,0x01, 0x18,0x03, 0x04,0x05, 0x06,0x07, 0x08,0x09, 0x0A,0x0B, 0x0C,0x0D, 0x0E,0x0F },  /* element 2 */
      { 0x00,0x01, 0x02,0x18, 0x04,0x05, 0x06,0x07, 0x08,0x09, 0x0A,0x0B, 0x0C,0x0D, 0x0E,0x0F },  /* element 3 */
      { 0x00,0x01, 0x02,0x03, 0x18,0x05, 0x06,0x07, 0x08,0x09, 0x0A,0x0B, 0x0C,0x0D, 0x0E,0x0F },  /* element 4 */
      { 0x00,0x01, 0x02,0x03, 0x04,0x18, 0x06,0x07, 0x08,0x09, 0x0A,0x0B, 0x0C,0x0D, 0x0E,0x0F },  /* element 5 */
      { 0x00,0x01, 0x02,0x03, 0x04,0x05, 0x18,0x07, 0x08,0x09, 0x0A,0x0B, 0x0C,0x0D, 0x0E,0x0F },  /* element 6 */
      { 0x00,0x01, 0x02,0x03, 0x04,0x05, 0x06,0x18, 0x08,0x09, 0x0A,0x0B, 0x0C,0x0D, 0x0E,0x0F },  /* element 7 */
      { 0x00,0x01, 0x02,0x03, 0x04,0x05, 0x06,0x07, 0x18,0x09, 0x0A,0x0B, 0x0C,0x0D, 0x0E,0x0F },  /* element 8 */
      { 0x00,0x01, 0x02,0x03, 0x04,0x05, 0x06,0x07, 0x08,0x18, 0x0A,0x0B, 0x0C,0x0D, 0x0E,0x0F },  /* element 9 */
      { 0x00,0x01, 0x02,0x03, 0x04,0x05, 0x06,0x07, 0x08,0x09, 0x18,0x0B, 0x0C,0x0D, 0x0E,0x0F },  /* element 10 */
      { 0x00,0x01, 0x02,0x03, 0x04,0x05, 0x06,0x07, 0x08,0x09, 0x0A,0x18, 0x0C,0x0D, 0x0E,0x0F },  /* element 11 */
      { 0x00,0x01, 0x02,0x03, 0x04,0x05, 0x06,0x07, 0x08,0x09, 0x0A,0x0B, 0x18,0x0D, 0x0E,0x0F },  /* element 12 */
      { 0x00,0x01, 0x02,0x03, 0x04,0x05, 0x06,0x07, 0x08,0x09, 0x0A,0x0B, 0x0C,0x18, 0x0E,0x0F },  /* element 13 */
      { 0x00,0x01, 0x02,0x03, 0x04,0x05, 0x06,0x07, 0x08,0x09, 0x0A,0x0B, 0x0C,0x0D, 0x18,0x0F },  /* element 14 */
      { 0x00,0x01, 0x02,0x03, 0x04,0x05, 0x06,0x07, 0x08,0x09, 0x0A,0x0B, 0x0C,0x0D, 0x0E,0x18 }   /* element 15 */
    #elif __BIG_ENDIAN__
      { 0x00,0x01, 0x02,0x03, 0x04,0x05, 0x06,0x07, 0x08,0x09, 0x0A,0x0B, 0x0C,0x0D, 0x0E,0x17 },  /* element 0 */
      { 0x00,0x01, 0x02,0x03, 0x04,0x05, 0x06,0x07, 0x08,0x09, 0x0A,0x0B, 0x0C,0x0D, 0x17,0x0F },  /* element 1 */
      { 0x00,0x01, 0x02,0x03, 0x04,0x05, 0x06,0x07, 0x08,0x09, 0x0A,0x0B, 0x0C,0x17, 0x0E,0x0F },  /* element 2 */
      { 0x00,0x01, 0x02,0x03, 0x04,0x05, 0x06,0x07, 0x08,0x09, 0x0A,0x0B, 0x17,0x0D, 0x0E,0x0F },  /* element 3 */
      { 0x00,0x01, 0x02,0x03, 0x04,0x05, 0x06,0x07, 0x08,0x09, 0x0A,0x17, 0x0C,0x0D, 0x0E,0x0F },  /* element 4 */
      { 0x00,0x01, 0x02,0x03, 0x04,0x05, 0x06,0x07, 0x08,0x09, 0x17,0x0B, 0x0C,0x0D, 0x0E,0x0F },  /* element 5 */
      { 0x00,0x01, 0x02,0x03, 0x04,0x05, 0x06,0x07, 0x08,0x17, 0x0A,0x0B, 0x0C,0x0D, 0x0E,0x0F },  /* element 6 */
      { 0x00,0x01, 0x02,0x03, 0x04,0x05, 0x06,0x07, 0x17,0x09, 0x0A,0x0B, 0x0C,0x0D, 0x0E,0x0F },  /* element 7 */
      { 0x00,0x01, 0x02,0x03, 0x04,0x05, 0x06,0x17, 0x08,0x09, 0x0A,0x0B, 0x0C,0x0D, 0x0E,0x0F },  /* element 8 */
      { 0x00,0x01, 0x02,0x03, 0x04,0x05, 0x17,0x07, 0x08,0x09, 0x0A,0x0B, 0x0C,0x0D, 0x0E,0x0F },  /* element 9 */
      { 0x00,0x01, 0x02,0x03, 0x04,0x17, 0x06,0x07, 0x08,0x09, 0x0A,0x0B, 0x0C,0x0D, 0x0E,0x0F },  /* element 10 */
      { 0x00,0x01, 0x02,0x03, 0x17,0x05, 0x06,0x07, 0x08,0x09, 0x0A,0x0B, 0x0C,0x0D, 0x0E,0x0F },  /* element 11 */
      { 0x00,0x01, 0x02,0x17, 0x04,0x05, 0x06,0x07, 0x08,0x09, 0x0A,0x0B, 0x0C,0x0D, 0x0E,0x0F },  /* element 12 */
      { 0x00,0x01, 0x17,0x03, 0x04,0x05, 0x06,0x07, 0x08,0x09, 0x0A,0x0B, 0x0C,0x0D, 0x0E,0x0F },  /* element 13 */
      { 0x00,0x17, 0x02,0x03, 0x04,0x05, 0x06,0x07, 0x08,0x09, 0x0A,0x0B, 0x0C,0x0D, 0x0E,0x0F },  /* element 14 */
      { 0x17,0x01, 0x02,0x03, 0x04,0x05, 0x06,0x07, 0x08,0x09, 0x0A,0x0B, 0x0C,0x0D, 0x0E,0x0F }   /* element 15 */
    #endif
  };
  __m128i_union t;
  t.as_long_long[__vsr_left_half_long_long_in_memory] = scalar;
  return (__m128i) vec_perm (v, t.as_m128i, permute_selector[element_from_right]);
}

/******************************************************* Extract ******************************************************/

/* Extract 32-bit int */
VECLIB_INLINE int vec_extract1swfrom4sw (__m128i from, const intlit2 element_from_right)
{
  __m128i_union from_union;
  from_union.as_m128i = from;
  #ifdef __LITTLE_ENDIAN__
    return from_union.as_int[element_from_right & 0x3];
  #elif __BIG_ENDIAN__
    return from_union.as_int[3 - (element_from_right & 0x3)];
  #endif
}

/* Extract 64-bit long long */
VECLIB_INLINE long long vec_extract1sdfrom2sd (__m128i from, const intlit1 element_from_right)
{
  __m128i_union from_union;
  from_union.as_m128i = from;
  #ifdef __LITTLE_ENDIAN__
    return from_union.as_long_long[element_from_right & 0x1];
  #elif __BIG_ENDIAN__
    return from_union.as_long_long[1 - (element_from_right & 0x1)];
  #endif
}

/* Extract 16-bit short from one of 8 16-bit shorts */
VECLIB_INLINE int vec_extract8sh (__m128i v, intlit3 element_from_right)
{
  __m128i_union t;
  #ifdef __LITTLE_ENDIAN__
    t.as_m128i = v;
    return t.as_short[element_from_right & 7];
  #elif __BIG_ENDIAN__
    static const vector unsigned char permute_selector[8] = {
    /* To extract specified halfword element into lowest halfword of the left half with other elements zeroed */
      { 0x00,0x01, 0x02,0x03, 0x04,0x05, 0x1E,0x1F, 0x08,0x09, 0x0A,0x0B, 0x0C,0x0D, 0x0E,0x0F },  /* element 0 */
      { 0x00,0x01, 0x02,0x03, 0x04,0x05, 0x1C,0x1D, 0x08,0x09, 0x0A,0x0B, 0x0C,0x0D, 0x0E,0x0F },  /* element 1 */
      { 0x00,0x01, 0x02,0x03, 0x04,0x05, 0x1A,0x1B, 0x08,0x09, 0x0A,0x0B, 0x0C,0x0D, 0x0E,0x0F },  /* element 2 */
      { 0x00,0x01, 0x02,0x03, 0x04,0x05, 0x18,0x19, 0x08,0x09, 0x0A,0x0B, 0x0C,0x0D, 0x0E,0x0F },  /* element 3 */
      { 0x00,0x01, 0x02,0x03, 0x04,0x05, 0x16,0x17, 0x08,0x09, 0x0A,0x0B, 0x0C,0x0D, 0x0E,0x0F },  /* element 4 */
      { 0x00,0x01, 0x02,0x03, 0x04,0x05, 0x14,0x15, 0x08,0x09, 0x0A,0x0B, 0x0C,0x0D, 0x0E,0x0F },  /* element 5 */
      { 0x00,0x01, 0x02,0x03, 0x04,0x05, 0x12,0x13, 0x08,0x09, 0x0A,0x0B, 0x0C,0x0D, 0x0E,0x0F },  /* element 6 */
      { 0x00,0x01, 0x02,0x03, 0x04,0x05, 0x10,0x11, 0x08,0x09, 0x0A,0x0B, 0x0C,0x0D, 0x0E,0x0F }   /* element 7 */
    };
    t.as_m128i = vec_perm (vec_splats ((unsigned char) 0), v, permute_selector[element_from_right & 7]);
    return (short) t.as_long_long[__vsr_left_half_long_long_in_memory];
  #endif
}

/* Extract upper bit of 16 8-bit chars */
VECLIB_INLINE int vec_extractupperbit16sb (__m128i v)
{
  __m128i_union t;
  t.as_m128i = v;
  int result = 0;
  #ifdef __LITTLE_ENDIAN__
    result |= (t.as_char[15] & 0x80) << (15-7);
    result |= (t.as_char[14] & 0x80) << (14-7);
    result |= (t.as_char[13] & 0x80) << (13-7);
    result |= (t.as_char[12] & 0x80) << (12-7);
    result |= (t.as_char[11] & 0x80) << (11-7);
    result |= (t.as_char[10] & 0x80) << (10-7);
    result |= (t.as_char[9]  & 0x80) <<  (9-7);
    result |= (t.as_char[8]  & 0x80) <<  (8-7);
    result |= (t.as_char[7]  & 0x80);
    result |= (t.as_char[6]  & 0x80) >>  (7-6);
    result |= (t.as_char[5]  & 0x80) >>  (7-5);
    result |= (t.as_char[4]  & 0x80) >>  (7-4);
    result |= (t.as_char[3]  & 0x80) >>  (7-3);
    result |= (t.as_char[2]  & 0x80) >>  (7-2);
    result |= (t.as_char[1]  & 0x80) >>  (7-1);
    result |= (t.as_char[0]  & 0x80) >>   7;
  #elif __BIG_ENDIAN__
    result |= (t.as_char[0]  & 0x80) << (15-7);
    result |= (t.as_char[1]  & 0x80) << (14-7);
    result |= (t.as_char[2]  & 0x80) << (13-7);
    result |= (t.as_char[3]  & 0x80) << (12-7);
    result |= (t.as_char[4]  & 0x80) << (11-7);
    result |= (t.as_char[5]  & 0x80) << (10-7);
    result |= (t.as_char[6]  & 0x80) <<  (9-7);
    result |= (t.as_char[7]  & 0x80) <<  (8-7);
    result |= (t.as_char[8]  & 0x80);
    result |= (t.as_char[9]  & 0x80) >>  (7-6);
    result |= (t.as_char[10] & 0x80) >>  (7-5);
    result |= (t.as_char[11] & 0x80) >>  (7-4);
    result |= (t.as_char[12] & 0x80) >>  (7-3);
    result |= (t.as_char[13] & 0x80) >>  (7-2);
    result |= (t.as_char[14] & 0x80) >>  (7-1);
    result |= (t.as_char[15] & 0x80) >>   7;
  #endif
  return result;
}

/* Extract upper bit of 2 64-bit doubles */
VECLIB_INLINE int vec_extractupperbit2dp (__m128d v)
{
  __m128d_union t;
  t.as_m128d = v;
  int result = 0;

  #ifdef __LITTLE_ENDIAN__
    result |= (t.as_int[3] & 0x80000000) >> 30;
    result |= (t.as_int[1] & 0x80000000) >> 31;
  #elif __BIG_ENDIAN__
    result |= (t.as_int[0] & 0x80000000) >> 30;
    result |= (t.as_int[2] & 0x80000000) >> 31;
  #endif
  return result;
}

/* Extract lower 32-bit int */
VECLIB_INLINE int vec_extractlowersw (__m128i from) {
  __m128i_union tempVal; tempVal.as_m128i = from;
  #ifdef __LITTLE_ENDIAN__
    return tempVal.as_int[0];
  #elif __BIG_ENDIAN__
    return tempVal.as_int[3];
  #endif
}

/************************************************ Convert integer to integer ******************************************/

/* Convert 8+8 16-bit shorts to 16 8-bit chars with signed saturation */
VECLIB_INLINE __m128i vec_packs8hto16sb (__m128i left, __m128i right)
{
  #ifdef __LITTLE_ENDIAN__
    return (__m128i) vec_packs ((vector signed short) left, (vector signed short) right);
  #elif __BIG_ENDIAN__
    return (__m128i) vec_packs ((vector signed short) right, (vector signed short) left);
  #endif
}

/* Convert 8+8 16-bit shorts to 16 8-bit chars with signed saturation - deprecated - use previous function */
VECLIB_INLINE __m128i vec_packs88hto16sb (__m128i left, __m128i right)
{
  return vec_packs8hto16sb (left, right);
}

/* Convert 4+4 32-bit ints to 8 16-bit shorts with signed saturation */
VECLIB_INLINE __m128i vec_packs4wto8sh (__m128i left, __m128i right)
{
  #ifdef __LITTLE_ENDIAN__
    return (__m128i) vec_packs ((vector signed int) left, (vector signed int) right);
  #elif __BIG_ENDIAN__
    return (__m128i) vec_packs ((vector signed int) right, (vector signed int) left);
  #endif
}

/* Convert 4+4 32-bit ints to 8 16-bit shorts with signed saturation - deprecated - use previous function */
VECLIB_INLINE __m128i vec_packs44wto8sh (__m128i left, __m128i right)
{
  return vec_packs4wto8sh (left, right);
}

/******************************************** Convert floating-point to integer ***************************************/

/* Convert 4 32-bit floats to 4 32-bit ints with truncation */
VECLIB_INLINE __m128i vec_converttruncating4spto4sw (__m128 a)
{
  __m128_all_union org_union, ret_union;
  org_union.as_m128 = a;
  ret_union.as_int[0] = (int) org_union.as_float[0];
  ret_union.as_int[1] = (int) org_union.as_float[1];
  ret_union.as_int[2] = (int) org_union.as_float[2];
  ret_union.as_int[3] = (int) org_union.as_float[3];
  return ret_union.as_m128i;
}

/* Convert 4 32-bit floats to 4 32-bit ints */
VECLIB_INLINE __m128i vec_convert4spto4sw (__m128 from)
{
  __m128_union from_union;
  from_union.as_m128 = from;
  __m128i_union result;
  /* Round to nearest integer */
  result.as_int[0] = from_union.as_float[0] > 0 ? (int) (from_union.as_float[0] + 0.5) : (int) (from_union.as_float[0] - 0.5);
  result.as_int[1] = from_union.as_float[1] > 0 ? (int) (from_union.as_float[1] + 0.5) : (int) (from_union.as_float[1] - 0.5);
  result.as_int[2] = from_union.as_float[2] > 0 ? (int) (from_union.as_float[2] + 0.5) : (int) (from_union.as_float[2] - 0.5);
  result.as_int[3] = from_union.as_float[3] > 0 ? (int) (from_union.as_float[3] + 0.5) : (int) (from_union.as_float[3] - 0.5);
  return result.as_m128i;
}

/* Convert 2 64-bit doubles to 2 32-bit ints with truncation */
#ifdef VECLIB_VSX
VECLIB_INLINE __m128i vec_Convert2dpto2sw (__m128d from) {
  #ifdef __ibmxl__
    vector signed int zeroes = {0x00000000, 0x00000000, 0x00000000, 0x00000000 };
    static const vector unsigned char permute_selector = {
    #ifdef __LITTLE_ENDIAN__
      0x10, 0x11, 0x12, 0x13, 0x1C, 0x1D, 0x1E, 0x1F, 0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07
    #elif __BIG_ENDIAN__
      0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x10, 0x11, 0x12, 0x13, 0x1C, 0x1D, 0x1E, 0x1F
    #endif
    };
    return (__m128i)vec_perm((vector signed int)zeroes , (vector signed int)vec_cts(from , 0) ,permute_selector );
  #else
    __m128i_union returnedVal;
    returnedVal.as_int[0] = (int)from[0];
    returnedVal.as_int[1] = (int)from[1];
    returnedVal.as_int[2] = 0;
    returnedVal.as_int[3] = 0;
    return returnedVal.as_m128i;
  #endif
}
#endif

/***************************************************** Arithmetic *****************************************************/

/* Add 16 8-bit chars */
VECLIB_INLINE __m128i vec_add16sb (__m128i left, __m128i right)
{
  return (__m128i) vec_add ((vector signed char) left, (vector signed char) right);
}

/* Add 8 16-bit shorts */
VECLIB_INLINE __m128i vec_add8sh (__m128i left, __m128i right)
{
  return (__m128i) vec_add ((vector signed short) left, (vector signed short) right);
}

/* Add 4 32-bit ints */
VECLIB_INLINE __m128i vec_add4sw (__m128i left, __m128i right)
{ return (__m128i) vec_add ((vector signed int) left, (vector signed int) right); }

/* Add 2 64-bit long longs */
VECLIB_INLINE __m128i vec_add2sd (__m128i left, __m128i right)
{
  return (__m128i) vec_add ((vector signed long long) left, (vector signed long long) right);
}

/* Add 16 8-bit chars with unsigned saturation */
VECLIB_INLINE __m128i vec_addsaturating16ub (__m128i left, __m128i right)
{ return (__m128i) vec_adds ((vector unsigned char) left, (vector unsigned char) right); }

/* Add 16 8-bit chars with signed saturation */
VECLIB_INLINE __m128i vec_addsaturating16sb (__m128i left, __m128i right)
{
  return (__m128i) vec_adds ((vector signed char) left, (vector signed char) right);
}

/* Add 8 16-bit shorts with signed saturation */
VECLIB_INLINE __m128i vec_addsaturating8sh (__m128i left, __m128i right)
{ return (__m128i) vec_adds ((vector signed short) left, (vector signed short) right); }

/* Add 8 16-bit shorts with unsigned saturation */
VECLIB_INLINE __m128i vec_addsaturating8uh (__m128i left, __m128i right)
{ return (__m128i) vec_adds ((vector unsigned short) left, (vector unsigned short) right); }

/* Subtract 16 8-bit chars */
VECLIB_INLINE __m128i vec_subtract16sb (__m128i left, __m128i right)
{
  return (__m128i) vec_sub ((vector signed char) left, (vector signed char) right);
}

/* Subtract 8 16-bit shorts */
VECLIB_INLINE __m128i vec_subtract8sh (__m128i left, __m128i right)
{
  return (__m128i) vec_sub ((vector signed short) left, (vector signed short) right);
}

/* Subtract 4 32-bit ints */
VECLIB_INLINE __m128i vec_subtract4sw (__m128i left, __m128i right)
{ return (__m128i) vec_sub ((vector signed int) left, (vector signed int) right); }

/* Subtract 2 64-bit long longs*/
VECLIB_INLINE __m128i vec_subtract2sd (__m128i left, __m128i right)
{
  return (__m128i) vec_sub ((vector signed long long) left, (vector signed long long) right);
}

/* Subtract 16 8-bit chars with unsigned saturation */
VECLIB_INLINE __m128i vec_subtractsaturating16ub (__m128i left, __m128i right)
{ return (__m128i) vec_subs ((vector unsigned char) left, (vector unsigned char) right); }

/* Subtract 16 8-bit chars with signed saturation */
VECLIB_INLINE __m128i vec_subtractsaturating16sb (__m128i left, __m128i right)
{
  return (__m128i) vec_subs ((vector signed char) left, (vector signed char) right);
}

/* Subtract 8 16-bit shorts with signed saturation */
VECLIB_INLINE __m128i vec_subtractsaturating8sh (__m128i left, __m128i right)
{ return (__m128i) vec_subs ((vector signed short) left, (vector signed short) right); }

/* Subtract 8 16-bit shorts with unsigned saturation */
VECLIB_INLINE __m128i vec_subtractsaturating8uh (__m128i left, __m128i right)
{ return (__m128i) vec_subs ((vector unsigned short) left, (vector unsigned short) right); }

/* Multiply lower 32-bit unsigned ints producing 2 64-bit unsigned long longs */
VECLIB_INLINE __m128i vec_multiplylower2uwto2ud (__m128i left, __m128i right)
{
  __m128i_union t;
  __m128i_union leftx;
  leftx.as_m128i = left;
  __m128i_union rightx;
  rightx.as_m128i = right;
  #ifdef __BIG_ENDIAN__
    t.as_long_long [1] = ((unsigned long long) leftx.as_unsigned_int[3]) * ((unsigned long long) rightx.as_unsigned_int[3]);
    t.as_long_long [0] = ((unsigned long long) leftx.as_unsigned_int[1]) * ((unsigned long long) rightx.as_unsigned_int[1]);
  #elif __LITTLE_ENDIAN__
    t.as_long_long [0] = ((unsigned long long) leftx.as_unsigned_int[0]) * ((unsigned long long) rightx.as_unsigned_int[0]);
    t.as_long_long [1] = ((unsigned long long) leftx.as_unsigned_int[2]) * ((unsigned long long) rightx.as_unsigned_int[2]);
  #endif
  return t.as_m128i;
}

/* Multiply 8 16-bit signed shorts */
VECLIB_INLINE __m128i vec_multiply8sh (__m128i left, __m128i right)
{
  __m128i zero = vec_splat4sw(0);
  return (__m128i) vec_mladd ((vector signed short)left, (vector signed short)right, (vector signed short)zero);
}

/* Average 16 8-bit unsigned chars */
VECLIB_INLINE __m128i vec_average16ub (__m128i left, __m128i right)
{
  return (__m128i) vec_avg ((vector unsigned char) left, (vector unsigned char) right);
}

/* Average 8 16-bit unsigned shorts */
VECLIB_INLINE __m128i vec_average8uh (__m128i left, __m128i right)
{
  return (__m128i) vec_avg ((vector unsigned short) left, (vector unsigned short) right);
}

/* Max 8 16-bit shorts */
VECLIB_INLINE __m128i vec_max8sh (__m128i left, __m128i right)
{ return (__m128i) vec_max ((vector signed short) left, (vector signed short) right); }

/* Max 16 8-bit unsigned chars */
VECLIB_INLINE __m128i vec_max16ub (__m128i left, __m128i right)
{ return (__m128i) vec_max ((vector unsigned char) left, (vector unsigned char) right); }

/* Min 16 8-bit unsigned chars */
VECLIB_INLINE __m128i vec_min16ub (__m128i left, __m128i right)
{ return (__m128i) vec_min ((vector unsigned char) left, (vector unsigned char) right); }

/* Min 8 16-bit shorts */
VECLIB_INLINE __m128i vec_min8sh (__m128i left, __m128i right)
{
  return (__m128i) vec_min ((vector signed short) left, (vector signed short) right);
}

/* Sum 2 octets of absolute differences of 16 8-bit unsigned chars into 2 64-bit long longs */
VECLIB_INLINE __m128i vec_sumabsdiffs16ub (__m128i left, __m128i right)
{
  vector unsigned char minimums = vec_min ((vector unsigned char) left, (vector unsigned char) right);
  vector unsigned char maximums = vec_max ((vector unsigned char) left, (vector unsigned char) right);
  vector unsigned char absolute_differences = vec_sub (maximums, minimums);
  vector unsigned int int_sums = vec_sum4s (absolute_differences, vec_splats (0u));
  union {
    vector unsigned int v;
    int val[4];
  } int_sumsU;
  int_sumsU.v =int_sums;
  #ifdef __ibmxl__
    #ifdef __LITTLE_ENDIAN__
      return (__m128i) (vector signed int) { int_sumsU.val[0]+int_sumsU.val[1], 0, int_sumsU.val[2]+int_sumsU.val[3], 0 };
    #elif __BIG_ENDIAN__
      return (__m128i) vec_sum2s ((vector signed int) int_sums, vec_splats (0));
    #endif
  #elif __GNUC__
    /* No non-OpenPOWER vec_sum2s */
    #ifdef __LITTLE_ENDIAN__
      return (__m128i) (vector signed int) { int_sumsU.val[0]+int_sumsU.val[1], 0, int_sumsU.val[2]+int_sumsU.val[3], 0 };
    #elif __BIG_ENDIAN__
      return (__m128i) (vector signed int) { 0, int_sumsU.val[0]+int_sumsU.val[1], 0, int_sumsU.val[2]+int_sumsU.val[3] };
    #endif
  #else
    #error Compiler not supported yet.
  #endif
}

/* Multiply 4 16-bit shorts then add adjacent pairs with saturation to 4 32-bit ints */
VECLIB_INLINE __m128i vec_summultiply4sh (__m128i left, __m128i right)
{
  /* Produce non-overflow product */
  vector signed int multiply_even = vec_mule ((vector signed short) left, (vector signed short) right);
  vector signed int multiply_odd = vec_mulo ((vector signed short) left, (vector signed short) right);
  /* Horizontal add */
  vector signed int sums = vec_add (multiply_even, multiply_odd);
  return (__m128i) sums;
}

/* Absolute value 16 8-bit chars */
VECLIB_INLINE __m128i vec_Abs16sb (__m128i a) {
  return (__m128i) vec_abs ((vector signed char) a);
}

/* Absolute value 8 16-bit shorts */
VECLIB_INLINE __m128i vec_Abs8sh (__m128i a) {
  return (__m128i) vec_abs ((vector signed short) a);
}

/* Absolute value 4 32-bit ints */
VECLIB_INLINE __m128i vec_Abs4sw (__m128i a) {
  return (__m128i) vec_abs((vector signed int)a);
}

#ifdef __LITTLE_ENDIAN__
  #define LEleft_BEright left
  #define LEright_BEleft right
#elif __BIG_ENDIAN__
  #define LEleft_BEright right
  #define LEright_BEleft left
#endif

/* Horizontally add 4+4 adjacent pairs of 16-bit shorts to 8 16-bit shorts - (a0+a1, a2+a3, a4+a5, a6+a7, b0+b1, b2+b3, b4+b5, b6+b7) */
VECLIB_INLINE __m128i vec_horizontalAdd8sh (__m128i left, __m128i right) {
  static vector unsigned char permuteOdds = (vector unsigned char) {
  #ifdef __LITTLE_ENDIAN__
    0x02,0x03,0x06,0x07, 0x0A,0x0B,0x0E,0x0F, 0x12,0x13,0x16,0x17, 0x1A,0x1B,0x1E,0x1F
  #elif __BIG_ENDIAN__
    0x00,0x01,0x04,0x05, 0x08,0x09,0x0C,0x0D, 0x10,0x11,0x14,0x15, 0x18,0x19,0x1C,0x1D
  #endif
  };
  return (__m128i) vec_add(
    (vector signed short)vec_perm(LEleft_BEright, LEright_BEleft, permuteOdds),
    (vector signed short)vec_pack((vector unsigned int)LEleft_BEright, (vector unsigned int)LEright_BEleft));
}

/* Horizontally add 2+2 adjacent pairs of 32-bit ints to 4 32-bit ints */
VECLIB_INLINE __m128i vec_partialhorizontaladd2sw (__m128i left, __m128i right)
{
  #ifdef __LITTLE_ENDIAN__
    static vector unsigned char addend_1_permute_mask = (vector unsigned char)
      { 0x00,0x01,0x02,0x03, 0x08,0x09,0x0A,0x0B, 0x10,0x11,0x12,0x13, 0x18,0x19,0x1A,0x1B };
    static vector unsigned char addend_2_permute_mask = (vector unsigned char)
      { 0x04,0x05,0x06,0x07, 0x0C,0x0D,0x0E,0x0F, 0x14,0x15,0x16,0x17, 0x1C,0x1D,0x1E,0x1F };
  #elif __BIG_ENDIAN__
    static vector unsigned char addend_1_permute_mask = (vector unsigned char)
      { 0x14,0x15,0x16,0x17, 0x1C,0x1D,0x1E,0x1F, 0x04,0x05,0x06,0x07, 0x0C,0x0D,0x0E,0x0F };
    static vector unsigned char addend_2_permute_mask = (vector unsigned char)
      { 0x10,0x11,0x12,0x13, 0x18,0x19,0x1A,0x1B, 0x00,0x01,0x02,0x03, 0x08,0x09,0x0A,0x0B };
  #endif
  vector signed int addend_1 = vec_perm ((vector signed int) left, (vector signed int) right, addend_1_permute_mask);
  vector signed int addend_2 = vec_perm ((vector signed int) left, (vector signed int) right, addend_2_permute_mask);
  return (__m128i) vec_add (addend_1, addend_2);
}

/* Horizontally add 4+4 adjacent pairs of 16-bit shorts to 8 16-bit shorts with saturation - (a0+a1, a2+a3, a4+a5, a6+a7, b0+b1, b2+b3, b4+b5, b6+b7) */
VECLIB_INLINE __m128i vec_horizontalAddsaturating8sh (__m128i left, __m128i right) {
  static vector unsigned char permuteOdds = (vector unsigned char) {
    #ifdef __LITTLE_ENDIAN__
      0x02,0x03,0x06,0x07, 0x0A,0x0B,0x0E,0x0F, 0x12,0x13,0x16,0x17, 0x1A,0x1B,0x1E,0x1F
    #elif __BIG_ENDIAN__
      0x00,0x01,0x04,0x05, 0x08,0x09,0x0C,0x0D, 0x10,0x11,0x14,0x15, 0x18,0x19,0x1C,0x1D
    #endif
  };
  return (__m128i) vec_adds(
    (vector signed short)vec_perm(LEleft_BEright, LEright_BEleft, permuteOdds),
    (vector signed short)vec_pack((vector unsigned int)LEleft_BEright, (vector unsigned int)LEright_BEleft));
}

/* Horizontally subtract 4+4 adjacent pairs of 16-bit shorts to 8 16-bit shorts - (a0-a1, a2-a3, a4-a5, a6-a7, b0-b1, b2-b3, b4-b5, b6-b7) */
VECLIB_INLINE __m128i vec_horizontalSub8sh (__m128i left, __m128i right) {
  static vector unsigned char permuteOdds = (vector unsigned char) {
    #ifdef __LITTLE_ENDIAN__
      0x02,0x03,0x06,0x07, 0x0A,0x0B,0x0E,0x0F, 0x12,0x13,0x16,0x17, 0x1A,0x1B,0x1E,0x1F
    #elif __BIG_ENDIAN__
      0x00,0x01,0x04,0x05, 0x08,0x09,0x0C,0x0D, 0x10,0x11,0x14,0x15, 0x18,0x19,0x1C,0x1D
    #endif
  };
  return (__m128i) vec_sub(
    (vector signed short)vec_pack((vector unsigned int)LEleft_BEright, (vector unsigned int)LEright_BEleft),
    (vector signed short)vec_perm(LEleft_BEright, LEright_BEleft, permuteOdds)
  );
}

/* Horizontally subtract 2+2 adjacent pairs of 32-bit ints to 4 32-bit ints */
VECLIB_INLINE __m128i vec_partialhorizontalsubtract2sw (__m128i left, __m128i right)
{
  #ifdef __LITTLE_ENDIAN__
    static vector unsigned char subend_1_permute_mask = (vector unsigned char)
      { 0x00,0x01,0x02,0x03, 0x08,0x09,0x0A,0x0B, 0x10,0x11,0x12,0x13, 0x18,0x19,0x1A,0x1B };
    static vector unsigned char subend_2_permute_mask = (vector unsigned char)
      { 0x04,0x05,0x06,0x07, 0x0C,0x0D,0x0E,0x0F, 0x14,0x15,0x16,0x17, 0x1C,0x1D,0x1E,0x1F };
  #elif __BIG_ENDIAN__
    static vector unsigned char subend_1_permute_mask = (vector unsigned char)
      { 0x14,0x15,0x16,0x17, 0x1C,0x1D,0x1E,0x1F, 0x04,0x05,0x06,0x07, 0x0C,0x0D,0x0E,0x0F };
    static vector unsigned char subend_2_permute_mask = (vector unsigned char)
      { 0x10,0x11,0x12,0x13, 0x18,0x19,0x1A,0x1B, 0x00,0x01,0x02,0x03, 0x08,0x09,0x0A,0x0B };
  #endif
  vector signed int subend_1 = vec_perm ((vector signed int) left, (vector signed int) right, subend_1_permute_mask);
  vector signed int subend_2 = vec_perm ((vector signed int) left, (vector signed int) right, subend_2_permute_mask);
  return (__m128i) vec_sub (subend_1, subend_2);
}

/* Horizontally subtract 4+4 adjacent pairs of 16-bit shorts to 8 16-bit shorts with saturation - (a0-a1, a2-a3, a4-a5, a6-a7, b0-b1, b2-b3, b4-b5, b6-b7) */
VECLIB_INLINE __m128i vec_horizontalSubtractsaturating8sh (__m128i left, __m128i right) {
  static vector unsigned char permuteOdds = (vector unsigned char) {
    #ifdef __LITTLE_ENDIAN__
      0x02,0x03,0x06,0x07, 0x0A,0x0B,0x0E,0x0F, 0x12,0x13,0x16,0x17, 0x1A,0x1B,0x1E,0x1F
    #elif __BIG_ENDIAN__
      0x00,0x01,0x04,0x05, 0x08,0x09,0x0C,0x0D, 0x10,0x11,0x14,0x15, 0x18,0x19,0x1C,0x1D
    #endif
  };
  return (__m128i) vec_subs(
    (vector signed short)vec_pack((vector unsigned int)LEleft_BEright, (vector unsigned int)LEright_BEleft),
    (vector signed short)vec_perm(LEleft_BEright, LEright_BEleft, permuteOdds)
  );
}

/* Multiply 16 8-bit u*s chars then add adjacent 16-bit products with signed saturation */
VECLIB_INLINE __m128i vec_Multiply16sbthenhorizontalAddsaturating8sh (__m128i left, __m128i right) {
  vector signed short zeroUpperHalfOfInts = { 0x00FF, 0x00FF, 0x00FF, 0x00FF, 0x00FF, 0x00FF, 0x00FF, 0x00FF };
  //This extends the "left" to be 16-bit unsigned integer version of iteself. The vec_and prevents the makes simulates that the value is unsigned
  __m128i_union unionLowerLeft;  unionLowerLeft.as_vector_signed_short =  vec_and( vec_unpackl( (vector signed char)left ), zeroUpperHalfOfInts );
  __m128i_union unionUpperLeft;  unionUpperLeft.as_vector_signed_short =  vec_and( vec_unpackh( (vector signed char)left ), zeroUpperHalfOfInts );
  //This extends the "Right" to 16-bit signed integer version of iteself. The unpack will extend the sign to the entire integer
  __m128i_union unionLowerRight; unionLowerRight.as_vector_signed_short = vec_unpackl( (vector signed char)right );
  __m128i_union unionUpperRight; unionUpperRight.as_vector_signed_short = vec_unpackh( (vector signed char)right );

  __m128i_union unionResultLower;
  __m128i_union unionResultUpper;

  unionResultLower.as_vector_signed_int = vec_add(
    vec_mulo(unionLowerLeft.as_vector_signed_short, unionLowerRight.as_vector_signed_short),
    vec_mule(unionLowerLeft.as_vector_signed_short, unionLowerRight.as_vector_signed_short)
  );

  unionResultUpper.as_vector_signed_int = vec_add(
    vec_mulo(unionUpperLeft.as_vector_signed_short, unionUpperRight.as_vector_signed_short),
    vec_mule(unionUpperLeft.as_vector_signed_short, unionUpperRight.as_vector_signed_short)
  );

  __m128i_union unionResultFinal;
  unionResultFinal.as_vector_signed_short = vec_packs( unionResultUpper.as_vector_signed_int, unionResultLower.as_vector_signed_int );
  return unionResultFinal.as_vector_unsigned_char;
}

/* Multiply 8 16-bit shorts, shift right 14, add 1 and shift right 1 to 8 16-bit shorts */
VECLIB_INLINE __m128i vec_Multiply8shExtractUpper (__m128i left, __m128i right) {

  __m128i_union odds; odds.as_vector_signed_int = vec_mulo((vector signed short)left, (vector signed short)right);
  __m128i_union evens; evens.as_vector_signed_int = vec_mule((vector signed short)left, (vector signed short)right);

  static vector unsigned int addVector = (vector unsigned int) { 0x00004000, 0x00004000, 0x00004000, 0x00004000 };
  static vector unsigned int shiftVector = (vector unsigned int) { 0x01010101, 0x01010101, 0x01010101, 0x01010101 };

  odds.as_vector_unsigned_int = vec_sll(vec_add(odds.as_vector_unsigned_int, addVector), shiftVector);
  evens.as_vector_unsigned_int = vec_sll(vec_add(evens.as_vector_unsigned_int, addVector), shiftVector);

  static vector unsigned char permuteMask = (vector unsigned char) {
    #ifdef __LITTLE_ENDIAN__
      0x02, 0x03, 0x12, 0x13,  0x06, 0x07, 0x16, 0x17,  0x0A, 0x0B, 0x1A, 0x1B,  0x0E, 0x0F, 0x1E, 0x1F
    #elif __BIG_ENDIAN__
      0x00, 0x01, 0x10, 0x11,  0x04, 0x05, 0x14, 0x15,  0x08, 0x09, 0x18, 0x19,  0x0C, 0x0D, 0x1C, 0x1D
    #endif
  };
  return (__m128i)vec_perm(evens.as_vector_signed_int, odds.as_vector_signed_int, permuteMask);
}

/* Negate 16 8-bit chars when mask is negative, zero when zero, else copy */
VECLIB_INLINE __m128i vec_conditionalNegate16sb (__m128i left, __m128i right) {
  vector signed char zeroes = {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00};
    return (__m128i)vec_sel(
    vec_sel(
      (vector signed char)vec_sub(zeroes, (vector signed char)left),
      (vector signed char)left,
      vec_cmpgt((vector signed char)right, zeroes)),
    zeroes,
    vec_cmpeq((vector signed char)right, zeroes));

}

/* Negate 8 16-bit shorts when mask is negative, zero when zero, else copy */
VECLIB_INLINE __m128i vec_conditionalNegate8sh (__m128i left, __m128i right) {
  vector signed short zeroes = {0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, };
  return (__m128i)vec_sel(
  vec_sel(
    (vector signed short)vec_sub(zeroes, (vector signed short)left),
    (vector signed short)left,
    vec_cmpgt((vector signed short)right, zeroes)),
  zeroes,
  vec_cmpeq((vector signed short)right, zeroes));
}

/* Negate 4 32-bit ints when mask is negative, zero when zero, else copy */
VECLIB_INLINE __m128i vec_conditionalNegate4sw (__m128i left, __m128i right) {
  vector signed int zeroes = {0x00000000, 0x00000000, 0x00000000, 0x00000000 };
    return (__m128i)vec_sel(
    vec_sel(
      (vector signed int)vec_sub(zeroes, (vector signed int)left),
      (vector signed int)left,
      vec_cmpgt((vector signed int)right, zeroes)),
    zeroes,
    vec_cmpeq((vector signed int)right, zeroes));
}

/* Multiply 4 32-bit signed ints */
VECLIB_INLINE __m128i vec_multiply4sw (__m128i left, __m128i right)
{
  #ifdef __ibmxl__
    return (__m128i) vec_mul ((vector signed int) left, (vector signed int) right);
  #elif __GNUC__
    __m128i_union left_union; left_union.as_m128i = left;
    __m128i_union right_union; right_union.as_m128i =right;
    __m128i_union ret_union;
    ret_union.as_int[0] = left_union.as_int[0] * right_union.as_int[0];
    ret_union.as_int[1] = left_union.as_int[1] * right_union.as_int[1];
    ret_union.as_int[2] = left_union.as_int[2] * right_union.as_int[2];
    ret_union.as_int[3] = left_union.as_int[3] * right_union.as_int[3];
    return ret_union.as_m128i;
  #endif
}

/* Max 4 32-bit signed ints */
VECLIB_INLINE __m128i vec_Max4sw (__m128i left, __m128i right) {
  return (__m128i)vec_max((vector signed int)left, (vector signed int)right);
}

/* Min 4 32-bit signed ints */
VECLIB_INLINE __m128i vec_Min4sw (__m128i left, __m128i right) {
  return (__m128i)vec_min((vector signed int)left, (vector signed int)right);
}


/******************************************************** Boolean *****************************************************/

/* Bitwise 128-bit and */
VECLIB_INLINE __m128i vec_bitand1q (__m128i left, __m128i right)
{ return (__m128i) vec_and ((vector unsigned char) left, (vector unsigned char) right); }

/* Bitwise 128-bit and not (reversed) */
VECLIB_INLINE __m128i vec_bitandnotleft1q (__m128i left, __m128i right)
{ return (__m128i) vec_andc ((vector unsigned char) right, (vector unsigned char) left); }

/* Bitwise 128-bit or */
VECLIB_INLINE __m128i vec_bitor1q (__m128i left, __m128i right)
{ return (__m128i) vec_or ((vector unsigned char) left, (vector unsigned char) right); }

/* Bitwise 128-bit xor */
VECLIB_INLINE __m128i vec_bitxor1q (__m128i left, __m128i right)
{ return (__m128i) vec_xor ((vector unsigned char) left, (vector unsigned char) right); }

/**************************************************** Unpack **********************************************************/

/* Unpack 8+8 8-bit chars from high halves and interleave */
VECLIB_INLINE __m128i vec_unpackhigh8sb (__m128i left, __m128i right)
{
  static const vector unsigned char permute_selector = {
    #ifdef __LITTLE_ENDIAN__
      0x08, 0x18, 0x09, 0x19,  0x0A, 0x1A, 0x0B, 0x1B,  0x0C, 0x1C, 0x0D, 0x1D,  0x0E, 0x1E, 0x0F, 0x1F
    #elif __BIG_ENDIAN__
      0x10, 0x00, 0x11, 0x01,  0x12, 0x02, 0x13, 0x03,  0x14, 0x04, 0x15, 0x05,  0x16, 0x06, 0x17, 0x07
    #endif
  };
  return vec_perm (left, right, permute_selector);
}

/* Unpack 8+8 8-bit chars from high halves and interleave - deprecated - use previous function */
VECLIB_INLINE __m128i vec_unpackhigh88sb (__m128i left, __m128i right)
{
  return vec_unpackhigh8sb (left, right);
}

/* Unpack 4+4 16-bit shorts from high halves and interleave */
VECLIB_INLINE __m128i vec_unpackhigh4sh (__m128i left, __m128i right)
{
  static const vector unsigned char permute_selector = {
    #ifdef __LITTLE_ENDIAN__
      0x17, 0x16, 0x07, 0x06,  0x15, 0x14, 0x05, 0x04,  0x13, 0x12, 0x03, 0x02,  0x11, 0x10, 0x01, 0x00
    #elif __BIG_ENDIAN__
      0x10, 0x11, 0x00, 0x01,  0x12, 0x13, 0x02, 0x03,  0x14, 0x15, 0x04, 0x05,  0x16, 0x17, 0x06, 0x07
    #endif
  };
  return vec_perm (left, right, permute_selector);
}

/* Unpack 4+4 16-bit shorts from high halves and interleave - deprecated - use previous function */
VECLIB_INLINE __m128i vec_unpackhigh44sh (__m128i left, __m128i right)
{
  return vec_unpackhigh4sh (left, right);
}

/* Unpack 8+8 8-bit chars from low halves and interleave */
VECLIB_INLINE __m128i vec_unpacklow8sb (__m128i left, __m128i right)
{
  static const vector unsigned char permute_selector = {
    #ifdef __LITTLE_ENDIAN__
      0x00, 0x10, 0x01, 0x11, 0x02, 0x12, 0x03, 0x13, 0x04, 0x14, 0x05, 0x15, 0x06, 0x16, 0x07, 0x17
    #elif __BIG_ENDIAN__
      0x18, 0x08, 0x19, 0x09,  0x1A, 0x0A, 0x1B, 0x0B,  0x1C, 0x0C, 0x1D, 0x0D,  0x1E, 0x0E, 0x1F, 0x0F
    #endif
  };
  return vec_perm (left, right, permute_selector);
}

/* Unpack 8+8 8-bit chars from low halves and interleave - deprecated - use previous function */
VECLIB_INLINE __m128i vec_unpacklow88sb (__m128i left, __m128i right)
{
  return vec_unpacklow8sb (left, right);
}

/* Unpack 4+4 16-bit shorts from low halves and interleave */
VECLIB_INLINE __m128i vec_unpacklow4sh (__m128i left, __m128i right)
{
  static const vector unsigned char permute_selector = {
    #ifdef __LITTLE_ENDIAN__
      0x1F, 0x1E, 0x0F, 0x0E,  0x1D, 0x1C, 0x0D, 0x0C,  0x1B, 0x1A, 0x0B, 0x0A,  0x19, 0x18, 0x09, 0x08
    #elif __BIG_ENDIAN__
      0x18, 0x19, 0x08, 0x09,  0x1A, 0x1B, 0x0A, 0x0B,  0x1C, 0x1D, 0x0C, 0x0D,  0x1E, 0x1F, 0x0E, 0x0F
    #endif
  };
  return vec_perm (left, right, permute_selector);
}

/* Unpack 4+4 16-bit shorts from low halves and interleave - deprecated - use previous function */
VECLIB_INLINE __m128i vec_unpacklow44sh (__m128i left, __m128i right)
{
  return vec_unpacklow4sh (left, right);
}

/* Unpack 2+2 32-bit ints from low halves and interleave */
VECLIB_INLINE __m128i vec_unpacklow2sw (__m128i to_even, __m128i to_odd)
{
  static const vector unsigned char permute_selector = {
    #ifdef __LITTLE_ENDIAN__
      0x00, 0x01, 0x02, 0x03,  0x10, 0x11, 0x12, 0x13,  0x04, 0x05, 0x06, 0x07,  0x14, 0x15, 0x16, 0x17
    #elif __BIG_ENDIAN__
      0x18, 0x19, 0x1A, 0x1B,  0x08, 0x09, 0x0A, 0x0B,  0x1C, 0x1D, 0x1E, 0x1F,  0x0C, 0x0D, 0x0E, 0x0F
    #endif
  };
  return vec_perm (to_even, to_odd, permute_selector);
}

/* Unpack 2+2 32-bit ints from high halves and interleave */
VECLIB_INLINE __m128i vec_unpackhigh2sw (__m128i to_even, __m128i to_odd)
{
  static const vector unsigned char permute_selector = {
    #ifdef __LITTLE_ENDIAN__
      0x08, 0x09, 0x0A, 0x0B,  0x18, 0x19, 0x1A, 0x1B,  0x0C, 0x0D, 0x0E, 0x0F,  0x1C, 0x1D, 0x1E, 0x1F
    #elif __BIG_ENDIAN__
      0x10, 0x11, 0x12, 0x13,  0x00, 0x01, 0x02, 0x03,  0x14, 0x15, 0x16, 0x17,  0x04, 0x05, 0x06, 0x07
    #endif
  };
  return vec_perm (to_even, to_odd, permute_selector);
}

/* Unpack 1+1 64-bit long longs from low halves and interleave */
VECLIB_INLINE __m128i vec_unpacklow1sd (__m128i to_even, __m128i to_odd)
{
  static const vector unsigned char permute_selector = {
    #ifdef __LITTLE_ENDIAN__
      0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07,  0x10, 0x11, 0x12, 0x13, 0x14, 0x15, 0x16, 0x17
    #elif __BIG_ENDIAN__
      0x18, 0x19, 0x1A, 0x1B, 0x1C, 0x1D, 0x1E, 0x1F,  0x08, 0x09, 0x0A, 0x0B, 0x0C, 0x0D, 0x0E, 0x0F
    #endif
  };
  return vec_perm (to_even, to_odd, permute_selector);
}

/* Unpack 1+1 64-bit long longs from high halves and interleave */
VECLIB_INLINE __m128i vec_unpackhigh1sd (__m128i to_even, __m128i to_odd)
{
  static const vector unsigned char permute_selector = {
    #ifdef __LITTLE_ENDIAN__
      0x08, 0x09, 0x0A, 0x0B, 0x0C, 0x0D, 0x0E, 0x0F,  0x18, 0x19, 0x1A, 0x1B, 0x1C, 0x1D, 0x1E, 0x1F
    #elif __BIG_ENDIAN__
      0x10, 0x11, 0x12, 0x13, 0x14, 0x15, 0x16, 0x17,  0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07
    #endif
  };
  return vec_perm (to_even, to_odd, permute_selector);
}

/****************************************************** Shift *********************************************************/

/*- SSE2 shifts >= 32 produce 0; Altivec/MMX shifts by count%count_size. */
/*- The Altivec spec says the element shift count vector register must have a shift count in each element */
/*- and the shift counts may be different for each element. */
/*- The PowerPC Architecture says all elements must contain the same shift count. That wins. */
/*- The element shift count_size is: byte shift: 3 bits (0-7), halfword: 4 bits (0-15), word: 5 bits (0-31). */
/*- For full vector shifts the Altivec/PowerPC bit shift count is in the rightmost 7 bits, */
/*- with a 4 bit slo/sro byte shift count then a 3 bit sll/srl bit shift count. */

/* Shift left */

/* Shift 8 16-bit shorts left logical */
VECLIB_INLINE __m128i vec_shiftleft8sh (__m128i v, __m128i count)
{
  __m128i_union t;
  t.as_m128i = count;
  unsigned short counts = t.as_short[0];
  if (counts >= 16)
  {
    return (__m128i) vec_splats (0);
  }
  else
  {
    return (__m128i) vec_sl ((vector signed short) v, (vector unsigned short) count);
  }
}

/* Shift 4 32-bit ints left logical */
VECLIB_INLINE __m128i vec_shiftleft4sw (__m128i v, __m128i count)
{
  __m128i_union t;
  t.as_m128i = count;
  unsigned int counts = t.as_unsigned_int[0];
  if (counts >= 32)
  {
    return (__m128i) vec_splats (0);
  }
  else
  {
    return (__m128i) vec_sl ((vector signed int) v, (vector unsigned int) count);
  }
}

/* Shift 2 64-bit long longs left logical */
VECLIB_INLINE __m128i vec_shiftleft2sd (__m128i v, __m128i count)
{
  __m128i_union t;
  t.as_m128i = count;
  unsigned long long counts = t.as_long_long[0];
  if (counts >= (long long) 64)
  {
    return (__m128i) vec_splats (0);
  }
  else
  {
    return (__m128i) vec_sl ((vector signed long long) v, (vector unsigned long long) count);
  }
}

/* Shift 8 16-bit shorts left logical immediate */
VECLIB_INLINE __m128i vec_shiftleftimmediate8sh (__m128i v, intlit8 count)
{
  if ((unsigned long) count >= 16)
  {
    /* SSE2 shifts >= element_size or < 0 produce 0; Altivec/MMX shifts by count%element_size. */
    return (__m128i) vec_splats (0);
  } else if (count == 0) {
    return v;
  } else {
    /* The PowerPC Architecture says all shift count fields must contain the same shift count. */
    __m128i_union replicated_count;
    replicated_count.as_vector_signed_short = vec_splats ((short) count);
    return (__m128i) vec_sl ((vector signed short) v, replicated_count.as_vector_unsigned_short);
  }
}

/* Shift 4 32-bit ints left logical immediate */
VECLIB_INLINE __m128i vec_shiftleftimmediate4sw (__m128i v, intlit8 count)
{
  if ((unsigned long) count >= 32)
  {
    /* SSE2 shifts >= element_size or < 0 produce 0; Altivec/MMX shifts by count%element_size. */
    return (__m128i) vec_splats (0);
  } else if (count == 0) {
    return v;
  } else {
    /* The PowerPC Architecture says all shift count fields must contain the same shift count. */
    __m128i_union replicated_count;
    replicated_count.as_vector_signed_int = vec_splats ((int) count);
    return (__m128i) vec_sl ((vector signed int) v, replicated_count.as_vector_unsigned_int);
  }
}

/* Shift 2 64-bit long longs left logical immediate */
VECLIB_INLINE __m128i vec_shiftleftimmediate2sd (__m128i v, intlit8 count)
{
  if ((unsigned long) count >= 64)
  {
    /* SSE2 shifts >= element_size or < 0 produce 0; Altivec/MMX shifts by count%element_size. */
    return (__m128i) vec_splats (0);
  } else if (count == 0) {
    return v;
  } else {
    /* The PowerPC Architecture says all shift count fields must contain the same shift count. */
    /* On Power7 vec_slo (vslo) does use just the documented bits 121:124. */
    /* On Power7 vec_sll (vsll) uses the lower 3 bits of each byte instead (legal). */
    __m128i_union replicated_count;
    replicated_count.as_vector_unsigned_char = vec_splats ((unsigned char) count);
    long long m = 0xFFFFFFFFFFFFFFFF << count;
    __m128i_union mask;
    mask.as_long_long[0] = m;
    mask.as_long_long[1] = m;
    return vec_and (vec_sll (vec_slo (v, replicated_count.as_m128i), replicated_count.as_m128i), mask.as_vector_unsigned_char);
  }
}

/* Shift 128-bits left logical immediate by bytes */
VECLIB_INLINE __m128i vec_shiftleftbytes1q (__m128i v, intlit8 bytecount)
{
  if ((unsigned long) bytecount >= 16)
  {
    /* SSE2 shifts >= element_size or < 0 produce 0; Altivec/MMX shifts by bytecount%element_size. */
    return (__m128i) vec_splats (0);
  } else if (bytecount == 0) {
    return v;
  } else {
    /* The PowerPC byte shift count must be multiplied by 8. */
    /* It need not but can be replicated, which handles both LE and BE shift count positioning. */
    __m128i_union replicated_count;
    replicated_count.as_m128i = vec_splat16sb (bytecount << 3);
    return (__m128i) vec_slo (v, replicated_count.as_m128i);
  }
}

/* Shift right */

/* Shift 2 64-bit long longs right logical immediate */
VECLIB_INLINE __m128i vec_shiftright2sd (__m128i v, __m128i count)
{
  __m128i_union t;
  t.as_m128i = count;
  unsigned long long counts = t.as_long_long[0];
  if (counts >= (long long) 64)
  {
    return (__m128i) vec_splats (0);
  }
  else
  {
    return (__m128i) vec_sr ((vector signed long long) v, (vector unsigned long long) count);
  }
}

/* Shift 8 16-bit shorts right logical immediate */
VECLIB_INLINE __m128i vec_shiftrightimmediate8sh (__m128i v, intlit8 count)
{
  if ((unsigned long) count >= 16) {
    /* SSE2 shifts >= element_size or < 0 produce 0; Altivec/MMX shifts by count%element_size. */
    return (__m128i) vec_splats (0);
  } else if (count == 0) {
    return v;
  } else {
    /* The PowerPC Architecture says all shift count fields must contain the same shift count. */
    __m128i_union replicated_count;
    replicated_count.as_vector_signed_short = vec_splats ((short) count);
    return (__m128i) vec_sr ((vector signed short) v, replicated_count.as_vector_unsigned_short);
  }
}

/* Shift 4 32-bit ints right logical immediate */
VECLIB_INLINE __m128i vec_shiftrightimmediate4sw (__m128i v, intlit8 count)
{
  if ((unsigned long) count >= 32)
  {
    /* SSE2 shifts >= element_size or < 0 produce 0; Altivec/MMX shifts by count%element_size. */
    return (__m128i) vec_splats (0);
  } else if (count == 0) {
    return v;
  } else {
    /* The PowerPC Architecture says all shift count fields must contain the same shift count. */
    __m128i_union replicated_count;
    replicated_count.as_vector_signed_int = vec_splats ((int) count);
    return (__m128i) vec_sr ((vector signed int) v, replicated_count.as_vector_unsigned_int);
  }
}

/* Shift 128-bits right logical immediate by bytes */
VECLIB_INLINE __m128i vec_shiftrightbytes1q (__m128i v, intlit8 bytecount)
{
  if ((unsigned long) bytecount >= 16)
  {
    /* SSE2 shifts >= element_size or < 0 produce 0; Altivec/MMX shifts by bytecount%element_size. */
    return (__m128i) vec_splats (0);
  } else if (bytecount == 0) {
    return v;
  } else {
    /* The PowerPC byte shift count must be multiplied by 8. */
    /* It need not but can be replicated, which handles both LE and BE shift count positioning. */
    __m128i_union replicated_count;
    replicated_count.as_m128i = vec_splat16sb (bytecount << 3);
    /* AT gcc v7.1 may miscompile vec_sro as vec_slo? */
    return (__m128i) vec_sro (v, replicated_count.as_m128i);
  }
}

/* Shift 4 32-bit ints right arithmetic */
VECLIB_INLINE __m128i vec_shiftrightarithmetic4wimmediate (__m128i v, intlit8 count)
{
  if ((unsigned long) count >= 32)
  {
    /* SSE2 shifts >= element_size or < 0 produce 0; Altivec/MMX shifts by count%element_size. */
    return (__m128i) vec_splats (0);
  } else if (count == 0) {
    return v;
  } else {
    /* The PowerPC Architecture says all shift count fields must contain the same shift count. */
    __m128i_union replicated_count;
    replicated_count.as_vector_signed_int = vec_splats ((int) count);
    return (__m128i) vec_sra ((vector signed int) v, replicated_count.as_vector_unsigned_int);
  }
}

/* Shift 2 64-bit long longs right logical immediate */
VECLIB_INLINE __m128i vec_shiftrightlogical2dimmediate (__m128i v, intlit8 count)
{
  if ((unsigned long) count >= 64)
  {
    /* SSE2 shifts >= element_size or < 0 produce 0; Altivec/MMX shifts by count%element_size. */
    return (__m128i) vec_splats (0);
  } else if (count == 0) {
    return v;
  } else {
    /* The PowerPC Architecture says all shift count fields must contain the same shift count. */
    /* On Power7 vec_slo (vslo) does use just the documented bits 121:124. */
    /* On Power7 vec_sll (vsll) uses the lower 3 bits of each byte instead (legal). */
    __m128i_union replicated_count;
    replicated_count.as_vector_unsigned_char = vec_splats ((unsigned char) count);
    long long m = 0xFFFFFFFFFFFFFFFF >> count;
    __m128i_union mask;
    mask.as_long_long[0] = m;
    mask.as_long_long[1] = m;
    return vec_and (vec_srl (vec_sro (v, replicated_count.as_m128i), replicated_count.as_m128i), mask.as_vector_unsigned_char);
  }
}

/* Shift 8 16-bit shorts right arithmetic */
VECLIB_INLINE __m128i vec_shiftrightarithmetic8himmediate (__m128i v, intlit8 count)
{
  if ((unsigned long) count >= 16)
  {
    /* SSE2 shifts >= element_size or < 0 produce 0; Altivec/MMX shifts by count%element_size. */
    return (__m128i) vec_cmplt ((vector signed short) v, vec_splats ((signed short) 0));
  } else if (count == 0) {
    return v;
  } else {
    __m128i_union replicated_count;
    replicated_count.as_vector_signed_short = vec_splats ((signed short) count);
    return (__m128i) vec_sra ((vector signed short) v, replicated_count.as_vector_unsigned_short);
  }
}

/* Shift 4 32-bit ints left logical */
VECLIB_INLINE __m128i vec_shiftrightarithmetic4sw (__m128i v, __m128i count)
{
  __m128i_union t;
  t.as_m128i = count;
  #ifdef __LITTLE_ENDIAN__
    unsigned int counts = t.as_int[0];
  #elif __BIG_ENDIAN__
    unsigned int counts = t.as_int[3];
  #endif
  if (counts >= 32)
  {
    return (__m128i) vec_cmplt ((vector signed int) v, vec_splats ((signed int) 0));
  }
  else
  {
    return (__m128i) vec_sra ((vector signed int) v, vec_splats (counts));
  }
}

/* Shift 128+128-bits right into 128-bits */
VECLIB_INLINE __m128i vec_shiftright2dqw (__m128i left, __m128i right, int const count) {
  static const vector unsigned char zeroes = {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00};
  static const vector unsigned char permute_selector[16] = {
    #ifdef __LITTLE_ENDIAN__
      {0x10,0x11,0x12,0x13,0x14,0x15,0x16,0x17,0x18,0x19,0x1A,0x1B,0x1C,0x1D,0x1E,0x1F},
      {0x11,0x12,0x13,0x14,0x15,0x16,0x17,0x18,0x19,0x1A,0x1B,0x1C,0x1D,0x1E,0x1F,0x00},
      {0x12,0x13,0x14,0x15,0x16,0x17,0x18,0x19,0x1A,0x1B,0x1C,0x1D,0x1E,0x1F,0x00,0x01},
      {0x13,0x14,0x15,0x16,0x17,0x18,0x19,0x1A,0x1B,0x1C,0x1D,0x1E,0x1F,0x00,0x01,0x02},
      {0x14,0x15,0x16,0x17,0x18,0x19,0x1A,0x1B,0x1C,0x1D,0x1E,0x1F,0x00,0x01,0x02,0x03},
      {0x15,0x16,0x17,0x18,0x19,0x1A,0x1B,0x1C,0x1D,0x1E,0x1F,0x00,0x01,0x02,0x03,0x04},
      {0x16,0x17,0x18,0x19,0x1A,0x1B,0x1C,0x1D,0x1E,0x1F,0x00,0x01,0x02,0x03,0x04,0x05},
      {0x17,0x18,0x19,0x1A,0x1B,0x1C,0x1D,0x1E,0x1F,0x00,0x01,0x02,0x03,0x04,0x05,0x06},
      {0x18,0x19,0x1A,0x1B,0x1C,0x1D,0x1E,0x1F,0x00,0x01,0x02,0x03,0x04,0x05,0x06,0x07},
      {0x19,0x1A,0x1B,0x1C,0x1D,0x1E,0x1F,0x00,0x01,0x02,0x03,0x04,0x05,0x06,0x07,0x08},
      {0x1A,0x1B,0x1C,0x1D,0x1E,0x1F,0x00,0x01,0x02,0x03,0x04,0x05,0x06,0x07,0x08,0x09},
      {0x1B,0x1C,0x1D,0x1E,0x1F,0x00,0x01,0x02,0x03,0x04,0x05,0x06,0x07,0x08,0x09,0x0A},
      {0x1C,0x1D,0x1E,0x1F,0x00,0x01,0x02,0x03,0x04,0x05,0x06,0x07,0x08,0x09,0x0A,0x0B},
      {0x1D,0x1E,0x1F,0x00,0x01,0x02,0x03,0x04,0x05,0x06,0x07,0x08,0x09,0x0A,0x0B,0x0C},
      {0x1E,0x1F,0x00,0x01,0x02,0x03,0x04,0x05,0x06,0x07,0x08,0x09,0x0A,0x0B,0x0C,0x0D},
      {0x1F,0x00,0x01,0x02,0x03,0x04,0x05,0x06,0x07,0x08,0x09,0x0A,0x0B,0x0C,0x0D,0x0E}
    #elif __BIG_ENDIAN__
      {0x10,0x11,0x12,0x13,0x14,0x15,0x16,0x17,0x18,0x19,0x1A,0x1B,0x1C,0x1D,0x1E,0x1F},
      {0x0F,0x10,0x11,0x12,0x13,0x14,0x15,0x16,0x17,0x18,0x19,0x1A,0x1B,0x1C,0x1D,0x1E},
      {0x0E,0x0F,0x10,0x11,0x12,0x13,0x14,0x15,0x16,0x17,0x18,0x19,0x1A,0x1B,0x1C,0x1D},
      {0x0D,0x0E,0x0F,0x10,0x11,0x12,0x13,0x14,0x15,0x16,0x17,0x18,0x19,0x1A,0x1B,0x1C},
      {0x0C,0x0D,0x0E,0x0F,0x10,0x11,0x12,0x13,0x14,0x15,0x16,0x17,0x18,0x19,0x1A,0x1B},
      {0x0B,0x0C,0x0D,0x0E,0x0F,0x10,0x11,0x12,0x13,0x14,0x15,0x16,0x17,0x18,0x19,0x1A},
      {0x0A,0x0B,0x0C,0x0D,0x0E,0x0F,0x10,0x11,0x12,0x13,0x14,0x15,0x16,0x17,0x18,0x19},
      {0x09,0x0A,0x0B,0x0C,0x0D,0x0E,0x0F,0x10,0x11,0x12,0x13,0x14,0x15,0x16,0x17,0x18},
      {0x08,0x09,0x0A,0x0B,0x0C,0x0D,0x0E,0x0F,0x10,0x11,0x12,0x13,0x14,0x15,0x16,0x17},
      {0x07,0x08,0x09,0x0A,0x0B,0x0C,0x0D,0x0E,0x0F,0x10,0x11,0x12,0x13,0x14,0x15,0x16},
      {0x06,0x07,0x08,0x09,0x0A,0x0B,0x0C,0x0D,0x0E,0x0F,0x10,0x11,0x12,0x13,0x14,0x15},
      {0x05,0x06,0x07,0x08,0x09,0x0A,0x0B,0x0C,0x0D,0x0E,0x0F,0x10,0x11,0x12,0x13,0x14},
      {0x04,0x05,0x06,0x07,0x08,0x09,0x0A,0x0B,0x0C,0x0D,0x0E,0x0F,0x10,0x11,0x12,0x13},
      {0x03,0x04,0x05,0x06,0x07,0x08,0x09,0x0A,0x0B,0x0C,0x0D,0x0E,0x0F,0x10,0x11,0x12},
      {0x02,0x03,0x04,0x05,0x06,0x07,0x08,0x09,0x0A,0x0B,0x0C,0x0D,0x0E,0x0F,0x10,0x11},
      {0x01,0x02,0x03,0x04,0x05,0x06,0x07,0x08,0x09,0x0A,0x0B,0x0C,0x0D,0x0E,0x0F,0x10}
    #endif
  };

  if (count < 16) {
    return vec_perm(left, right, permute_selector[count]);
  }
  else if(count < 32) {
    return vec_perm(zeroes, left, permute_selector[count-16]);
  }
  else{
    return zeroes;
  }
}

/******************************************************* Permute ******************************************************/

/* Shuffle lower 4 16-bit shorts using mask, leaving upper half unchanged */
VECLIB_INLINE __m128i vec_permutelower4sh (__m128i v, intlit8 element_selectors)
{
  unsigned long element_selector_10 =  element_selectors       & 0x03;
  unsigned long element_selector_32 = (element_selectors >> 2) & 0x03;
  unsigned long element_selector_54 = (element_selectors >> 4) & 0x03;
  unsigned long element_selector_76 = (element_selectors >> 6) & 0x03;
  static const unsigned short permute_selectors [4] = {
    #ifdef __LITTLE_ENDIAN__
      0x0100, 0x0302, 0x0504, 0x0706
    #elif __BIG_ENDIAN__
      0x0E0F, 0x0C0D, 0x0A0B, 0x0809
    #endif
  };
  __m128i_union t;
  #ifdef __LITTLE_ENDIAN__
    t.as_long_long[__vsr_left_half_long_long_in_memory] = 0x0F0E0D0C0B0A0908;
    t.as_short[0] = permute_selectors[element_selector_10];
    t.as_short[1] = permute_selectors[element_selector_32];
    t.as_short[2] = permute_selectors[element_selector_54];
    t.as_short[3] = permute_selectors[element_selector_76];
  #elif __BIG_ENDIAN__
    t.as_long_long[__vsr_left_half_long_long_in_memory] = 0x0001020304050607;
    t.as_short[7] = permute_selectors[element_selector_10];
    t.as_short[6] = permute_selectors[element_selector_32];
    t.as_short[5] = permute_selectors[element_selector_54];
    t.as_short[4] = permute_selectors[element_selector_76];
  #endif
  return vec_perm (v, v, t.as_m128i);
}

/* Shuffle 4 32-bit ints using mask */
VECLIB_INLINE __m128i vec_permute4sw (__m128i v, intlit8 element_selectors)
{
  unsigned long element_selector_10 =  element_selectors       & 0x03;
  unsigned long element_selector_32 = (element_selectors >> 2) & 0x03;
  unsigned long element_selector_54 = (element_selectors >> 4) & 0x03;
  unsigned long element_selector_76 = (element_selectors >> 6) & 0x03;
  static const unsigned int permute_selectors [4] = {
    #ifdef __LITTLE_ENDIAN__
      0x03020100, 0x07060504, 0x0B0A0908, 0x0F0E0D0C
    #elif __BIG_ENDIAN__
      0x0C0D0E0F, 0x08090A0B, 0x04050607, 0x00010203
    #endif
  };
  __m128i_union t;
  #ifdef __LITTLE_ENDIAN__
    t.as_int[0] = permute_selectors[element_selector_10];
    t.as_int[1] = permute_selectors[element_selector_32];
    t.as_int[2] = permute_selectors[element_selector_54];
    t.as_int[3] = permute_selectors[element_selector_76];
  #elif __BIG_ENDIAN__
    t.as_int[3] = permute_selectors[element_selector_10];
    t.as_int[2] = permute_selectors[element_selector_32];
    t.as_int[1] = permute_selectors[element_selector_54];
    t.as_int[0] = permute_selectors[element_selector_76];
  #endif
  return vec_perm (v, v, t.as_m128i);
}

/* Shuffle 16 8-bit chars using mask */
VECLIB_INLINE __m128i vec_permute16sb (__m128i v, __m128i mask)
{
  /* Prevent out of range permute index */
  vector bool char zero_mask = vec_cmplt (mask, vec_splats ((unsigned char) 0x80));
  #ifdef __LITTLE_ENDIAN__
    vector unsigned char permute_vector = vec_and (mask, vec_splats ((unsigned char) 0xF));
  #elif __BIG_ENDIAN__
    vector unsigned char permute_vector = vec_and (mask, vec_splats ((unsigned char) 0xF));
    permute_vector = vec_sub (vec_splats ((unsigned char) 0xF), permute_vector);
  #endif
  vector unsigned char permute_result = vec_perm (v, v, permute_vector);
  return vec_and (permute_result, zero_mask);
}

/******************************************************* Compare ******************************************************/

/* Compare eq */

/* Compare 16 8-bit chars for == to vector mask */
VECLIB_INLINE __m128i vec_compareeq16sb (__m128i left, __m128i right)
{ return (__m128i) vec_cmpeq ((vector signed char) left, (vector signed char) right); }

/* Compare 8 16-bit shorts for == to vector mask */
VECLIB_INLINE __m128i vec_compareeq8sh (__m128i left, __m128i right)
{ return (__m128i) vec_cmpeq ((vector signed short) left, (vector signed short) right); }

/* Compare 4 32-bit ints for == to vector mask */
VECLIB_INLINE __m128i vec_compare4sw (__m128i left, __m128i right)
{ return (__m128i) vec_cmpeq ((vector signed int) left, (vector signed int) right); }

/* Compare lt */

/* Compare 16 8-bit chars for < to mask */
VECLIB_INLINE __m128i vec_comparelt16sb (__m128i left, __m128i right)
{
  return (__m128i) vec_cmplt ((vector signed char) left, (vector signed char) right);
}

/* Compare 8 16-bit shorts for < to vector mask */
VECLIB_INLINE __m128i vec_comparelt8sh (__m128i left, __m128i right)
{ return (__m128i) vec_cmplt ((vector signed short) left, (vector signed short) right); }

/* Compare 4 32-bit ints for < to vector mask */
VECLIB_INLINE __m128i vec_comparelt4sw (__m128i left, __m128i right)
{ return (__m128i) vec_cmplt ((vector signed int) left, (vector signed int) right); }

/* Compare gt */

/* Compare 16 8-bit signed chars for > to vector mask */
VECLIB_INLINE __m128i vec_comparegt16sb (__m128i left, __m128i right)
{
  return (__m128i) vec_cmpgt ((vector signed char) left, (vector signed char) right);
}

/* Compare 8 16-bit shorts for > to vector mask */
VECLIB_INLINE __m128i vec_comparegt8sh (__m128i left, __m128i right)
{ return (__m128i) vec_cmpgt ((vector signed short) left, (vector signed short) right); }

/* Compare 4 32-bit ints for > to vector mask */
VECLIB_INLINE __m128i vec_comparegt4sw (__m128i left, __m128i right)
{ return (__m128i) vec_cmpgt ((vector signed int) left, (vector signed int) right); }

/****************************************************** Cast **********************************************************/

/* Cast __m128 to __m128i */
VECLIB_INLINE __m128i vec_cast4spto1q (__m128 v)
{
  __m128_all_union v_union;
  v_union.as_m128 = v;
  return (__m128i) v_union.as_m128i;
}

/* Cast __m128d to __m128i */
VECLIB_INLINE __m128i vec_Cast2dpto4sw (__m128d from) {
  __m128_all_union newFrom; newFrom.as_m128d = from;
  return newFrom.as_m128i;
}
#endif
