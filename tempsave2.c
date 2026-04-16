



/* NOT WORKING WITH THE 32 BIT GEN POLYNOMIAL, THE OLD FILE IS SAVED IN TEMPSAVE.C WITH THE 28 BIT GEN POLYNOMIAL */






// #ifndef BCH_H
// #define BCH_H

// #include <stdio.h>
// #include <stdint.h>
// #include <stdbool.h>

// /* BCH(63, 32) Parameters */
// #define BCH_N 63      
// #define BCH_K 32      
// #define BCH_T 5       
// #define BCH_PARITY  31
// #define GF_PRIM_POLY  0x43



// void bch_init(void);

// uint64_t bch_encode(uint32_t data);

// void bch_compute_syndromes(uint64_t codeword, uint8_t syndromes[2 * BCH_T]);

// void bch_syn_error_locator(const uint8_t syndromes[2 * BCH_T], uint8_t sigma[BCH_T + 1]);

// uint64_t bch_syn_correction(uint64_t codeword, const uint8_t sigma[BCH_T + 1]);

// uint32_t bch_original_data(uint64_t corrected_codeword);

// uint32_t bch_decode(uint64_t codeword);





// #endif













#include "BCH.h"
#include <stdint.h>
#include <string.h>


/* 
 * GF(2^6) Internal Tables
 */

// @note: the lookup tables retuen the "i" i.e. the power in alpha^i 

static uint8_t gf_exp[126];  // antilog table {size == double of (2^6 - 1) to avoid mis matiching }
/* takes in alpha^i and returns the mapped binary value of the polynomial field element */

static uint8_t gf_log[64];   // log table
/* takes in the binary value of the polynomial field element to give out the corresponding mapped aplha^i value */

static uint64_t gen_poly;  


static inline uint8_t gf_pow(int i ) // what power of aplha does the element number(decimal) "n" map too?
{
    /* takes any decimal and give alpha^i for it from the lookup table */
    return gf_log[i];
}


static inline uint8_t gf_deci(int n) // what decimal element does alpha^n map too?
{
    /* takes in the "i" from alpha^i and returns its decimal equivalent */
    return gf_exp[n % 63];
}

static uint8_t gf_mul(uint8_t a, uint8_t b)
{
    /* MULTIPLIES two polynomial in binary, by converting them into their alpha^i value and adding, finally giving out the polynomial by reconverting */
    if(a == 0 || b == 0)
        return 0;

    return gf_exp[(gf_log[a] + gf_log[b]) % 63];
}

static uint64_t mod2_divide(uint64_t dividend, uint64_t divisor)
{
  
    // Find the bit length of divisor and dividend
    uint64_t div_len = 0, dis_len = 0;
    uint64_t temp;

    temp = divisor;

    while (temp) // finds the bit length
    {  
        div_len++;  
        temp >>= 1; 
    }

    /* VALUE OF GEN-POLYNOMIAL IS PRE CALCULATED TO BE 28 BITS SO NO NEED TO FIND AGAIN */
   

    temp = dividend;

    while (temp) // finds the bit length
    { 
        dis_len++; 
        temp >>= 1; 
    }

    /* VALUE OF GEN-POLYNOMIAL IS PRE CALCULATED TO BE 28 BITS SO NO NEED TO FIND AGAIN */

    // Perform modulo-2 division using XOR
    int64_t shift = dis_len - div_len;

    while (shift >= 0) {
        dividend ^= (divisor << shift);

        // Recalculate the length of the new dividend
        temp = dividend;
        dis_len = 0;

        while (temp) 
        { 
            dis_len++; 
            temp >>= 1; 
        }

        shift = dis_len - div_len;
    }

    return dividend; // The remainder
}



void bch_init(void)
{
    /* intializing function to create the GF(64) field elements */
    int i; 
    uint8_t x = 1;
    for (i = 0; i < 63; i++) {
        gf_exp[i] = x;
        gf_exp[i + 63] = x;
        gf_log[x] = (uint8_t)i;
        x <<= 1;
        if (x & 0x40){ 
            x ^= GF_PRIM_POLY;
        }
        x &= 0x3F; // masks the final value so that it stays in bound in 6-bit
    }
    gen_poly = 0x8018D62B3ULL; // pre-calculated

}

uint64_t bch_encode(uint32_t data)
{
    uint64_t msg = (uint64_t)data << BCH_PARITY;
    /* left shift the data by 31 bits to make space for the parity */
    /* but our gen-polynomial is a 28 bit i.e. degree 28 and dividing by it gives a remainder of 27 bits, thus left shift by 27 bits */
    /* this changes the bch algorithm from (63,32,5) to (59,32,5)*/
    /* note: in testing, the results were satisfactory with this algorithm */

    uint64_t parity = mod2_divide(msg, gen_poly) & ((1ULL << BCH_PARITY) - 1);
    /* parity is calculated using modulo-2 division of the shifted msg with the generator polynomial */
    /* note: in modulo-2 division, the subtraction is done using ex-or operation */
    return msg | parity;
}

void bch_compute_syndromes(uint64_t codeword, uint8_t syndromes[2 * BCH_T])
{
    int i, j;
    for (i = 1; i <= 2 * BCH_T; i++) {
        uint8_t s = 0;
        for (j = 0; j < BCH_N; j++)
            if ((codeword >> j) & 1)
                s ^= gf_deci(i * j);
        syndromes[i - 1] = s;
    }
}

void bch_syn_error_locator(const uint8_t syndromes[2 * BCH_T], uint8_t sigma[BCH_T + 1])
{
    uint8_t prev[BCH_T + 1], temp[BCH_T + 1];
    int L = 0, m = 1, i, j;
    uint8_t b = 1, d;

    memset(sigma, 0, (BCH_T + 1) * sizeof(uint8_t));
    memset(prev,  0, sizeof(prev));
    sigma[0] = prev[0] = 1;

    for (i = 0; i < 2 * BCH_T; i++) {
        d = syndromes[i];
        for (j = 1; j <= L; j++)
            d ^= gf_mul(sigma[j], syndromes[i - j]);

        if (d == 0) {
            m++;
        } else if (2 * L <= i) {
            memcpy(temp, sigma, sizeof(temp));
            for (j = m; j <= BCH_T; j++)
                sigma[j] ^= gf_mul(d, gf_mul(gf_deci(63 - gf_log[b]), prev[j - m]));
            L = i + 1 - L;
            memcpy(prev, temp, sizeof(prev));
            b = d; m = 1;
        } else {
            for (j = m; j <= BCH_T; j++)
                sigma[j] ^= gf_mul(d, gf_mul(gf_deci(63 - gf_log[b]), prev[j - m]));
            m++;
        }
    }
}

uint64_t bch_syn_correction(uint64_t codeword, const uint8_t sigma[BCH_T + 1])
{
    int degree = 0, i, j;
    for (i = BCH_T; i >= 1; i--)
        if (sigma[i]) { degree = i; break; }

    for (i = 0; i < BCH_N; i++) {
        uint8_t val = sigma[0];
        for (j = 1; j <= degree; j++)
            val ^= gf_mul(sigma[j], gf_deci(j * i));
        if (val == 0)
            codeword ^= (1ULL << ((BCH_N - i) % BCH_N));
    }
    return codeword;
}

uint32_t bch_original_data(uint64_t corrected_codeword)
{
    return (uint32_t)(corrected_codeword >> BCH_PARITY);
}


// uint32_t bch_decode(uint64_t codeword)
// {
//     uint64_t syndromes = bch_compute_syndromes(codeword);
//     uint64_t result = bch_syn_error_locator(syndromes);
//     uint64_t raw = bch_syn_correction(codeword, result);
//     return (uint32_t)(raw >> 31);
// }