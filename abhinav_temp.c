#include <stdio.h>
#include <stdint.h>
#include "BCH.h"

static uint8_t gf_exp[126];   /* Antilog table: alpha^i          */
static uint8_t gf_log[64];    /* Log table:     log_alpha(x)     */
static uint64_t gen_poly; 


// Modulo-2 division (XOR division)
// dividend and divisor are represented as integers (bit patterns)
// Returns the remainder
uint64_t mod2_divide(unsigned int dividend, unsigned int divisor)
{
    if (divisor == 0) {
        printf("Error: Divisor cannot be zero.\n");
        return 0;
    }

    // Find the bit length of divisor and dividend
    int div_len = 0, dis_len = 0;
    unsigned int temp;

    temp = divisor;
    while (temp) { div_len++;  temp >>= 1; } // finds the bit length

    temp = dividend;
    while (temp) { dis_len++; temp >>= 1; }// finds the bit length

    // If dividend is smaller than divisor, remainder is the dividend itself
    if (dis_len < div_len)
        return dividend;

    // Perform modulo-2 division using XOR
    int shift = dis_len - div_len;

    while (shift >= 0) {
        dividend ^= (divisor << shift);

        // Recalculate the length of the new dividend
        temp = dividend;
        dis_len = 0;
        while (temp) { dis_len++; temp >>= 1; }

        shift = dis_len - div_len;
    }

    return dividend; // The remainder
}

unsigned int mod2_divide_ans(unsigned int dividend, unsigned int divisor)
{
    if (divisor == 0) {
        printf("Error: Divisor cannot be zero.\n");
        return 0;
    }

    // Find bit lengths
    int div_len = 0, dis_len = 0;
    unsigned int temp;

    temp = divisor;
    while (temp) { div_len++; temp >>= 1; }

    temp = dividend;
    while (temp) { dis_len++; temp >>= 1; }

    // If dividend is smaller than divisor, quotient is 0
    if (dis_len < div_len)
        return 0;

    unsigned int quotient = 0;
    int shift = dis_len - div_len;

    while (shift >= 0) {
        // Set the corresponding bit in the quotient
        quotient |= (1 << shift);

        // XOR the dividend with the aligned divisor
        // wont have to worry about haven "n" zeros after XORing as we will recalculate the shift
        // the len will consider the msb the first "1" we get soo.
        dividend ^= (divisor << shift);

        // Recalculate bit length of new dividend
        temp = dividend;
        dis_len = 0;
        while (temp) { dis_len++; temp >>= 1; }

        shift = dis_len - div_len;
    }

    return quotient;
}

void init(void)
 {   int i;
    uint8_t x = 1;

    /* Build exp and log tables using primitive poly 0x43 (x^6 + x + 1) */
    for (i = 0; i < 63; i++)
    {
        gf_exp[i] = x;
        gf_exp[i + 63] = x;
        gf_log[x] = i;
        x <<= 1;
        if (x & 0x40) // check overflow
            x ^= GF_PRIM_POLY; // modulo-2 division to bring the value back in GF(2^6)
    }

    printf("primitive polynomial = %#x\n" , GF_PRIM_POLY);
    for (int i = 0; i<63; i++ ){
        printf("gf_exp[%d] == %d && gf_log[%d] == %d\n" , i , gf_exp[i] , i , gf_log[i]);
    }

}


/* =================================== TESTING STUFF ============================================================*/

static inline uint8_t gf_pow(int n )  // what power of aplha does the element number(decimal) "n" map too? 
{
    /* takes any decimal and give alpha^i for it from the lookup table */
    return gf_log[n];
}

static inline uint8_t gf_deci(int n) // what decimal element does alpha^n map too?
{
    /* takes in the "i" from alpha^i and returns its decimal equivalent */
    return gf_exp[n % 63];
}

static uint8_t gf_mul(uint8_t a, uint8_t b)
{
    if(a == 0 || b == 0)
        return 0;

    return gf_exp[(gf_log[a] + gf_log[b]) % 63];
}

static void test_generator_polynomial(void)
{
    /* for a gf(8), we use the primitive polynomial x^3 + x + 1 */
    /*
        important note: we use mod with the primitive polynomial as a wrap around condition, to prevent overflow from the GF(2^3) i.e. only 3 bits.
        by observation, we can see that mod of two polynomials in just EXoring in binary
        for example: let alpha be the primitive element and p(x) == x^3 + x + 1  (in binary GF(2) == 1011)

        PRIM_ELEMMENT      BINARY(2^3)                                                 POLYINOMIAL
 
        alpha ^ 0   -->        001                                                         1

        alpha ^ 1   -->        010                                                         x

        alpha ^ 2   -->        100                                                        (x)^2

        alpha ^ 3   -->      1 000 (overflown)
                            thus we operate with == (mod p(x))                            (x)^3 (mod p(x)) ==>
                            ==> 1000 ^ 1011 ==> 0011                                      x + 1 (exact remainder found using polynomial division)
                            and 0011 in polynomial is x + 1 

                               011                                                        x + 1

        alpha ^ 4   -->        110                                                       (x)^2 + 1

        alpha ^ 5   -->        111                                                       (x)^2 + x + 1

        alpha ^ 6   -->        101                                                        (x)^2 + 1

        alpha ^ 7   -->        001    == alpha ^ 0                                            1
                            ( cyclic manner now)
    */

    uint8_t PRIM_POLY_GF_8 = 0xb;
    uint8_t gf_log_8[8]; // takes in polynomial field elemnet in binary to give i from (alpha^i)
    uint8_t gf_antilog_8[8]; // takes in i from (alpha^i) and gives field element polynomial in binary

    uint8_t x = 1;

    for(int i = 0; i<8; i++)
    {
        gf_antilog_8[i] = x;
        gf_log[x] = i;

        x <<= 1;
        if(x & 0x8)
        {
            x ^= PRIM_POLY_GF_8;
        }

    }

     printf("primitive polynomial for gf(8) = %#x\n" , PRIM_POLY_GF_8);
    for (int i = 0; i<8; i++ ){
        printf("gf_antilog_8[%d] = %d\n" , i , gf_antilog_8[i]);
    }

    /* WORKS BBG*/
}



int main()
{
    
    unsigned int dividend = 0b11010011101100; // example binary pattern
    unsigned int divisor  = 0b1011;

    uint64_t remainder = mod2_divide(dividend, divisor);
    uint64_t ans = mod2_divide_ans(dividend, divisor);

    printf("Remainder: %#x (binary representation)\n", remainder);
    printf("ans: %#x (binary representation)\n", ans);

    init();

    int n = 66; // in GF(2^6), 66 === 3
    int a = 15;
    int b = 3;

    uint8_t x = gf_pow(n); // returns gf_exp[n] but with wrap around, so n is not limitted to 0-63
    /* as wrap around condition is used, (%63), all integeres are wraped
       gf_log[3] == 8, or in polynomial, alpha^8
       and gf_log[66] == 8,
    */

   
   printf("gf_pow(%d) = %d\n" , a ,gf_pow(a));
   printf("gf_pow(%d) = %d\n" , b ,gf_pow(b));
   printf("gf_deci(%d) = %d\n" , gf_pow(a) + gf_pow(b), gf_deci(gf_pow(a) + gf_pow(b)));
   printf("gf_mul(%d,%d) = %d\n" ,a,b, gf_mul(a,b));

   test_generator_polynomial();

   unsigned long long l = 0x86e8113ULL;
    int bits = 0;

    while (l) {
        bits++;
        l >>= 1;
    }

    printf("bit length = %d\n", bits);
   
    return 0;
}