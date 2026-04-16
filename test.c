#include "BCH.h"

#define DATA 0x126ff5ea

int main()
{
    bch_init();

    uint64_t codeword = bch_encode(DATA);
    printf("original data = %#x\n",   DATA);
    printf("codeword = %#llx\n", (unsigned long long)codeword);

    uint64_t error_codeword = codeword;
    error_codeword ^= (1ULL << 10);
    error_codeword ^= (1ULL << 5);
    error_codeword ^= (1ULL << 32);
    error_codeword ^= (1ULL << 18);
    error_codeword ^= (1ULL << 60);
    error_codeword ^= (1ULL << 43);
    error_codeword ^= (1ULL << 40);
    printf("error codeword = %#llx\n", (unsigned long long)error_codeword);

    uint8_t syn[2 * BCH_T];
    bch_compute_syndromes(error_codeword, syn);
    printf("syndromes = ");
    for (int i = 0; i < 2*BCH_T; i++) printf("%02x ", syn[i]);
    printf("\n");

    uint8_t sigma[BCH_T + 1];
    bch_syn_error_locator(syn, sigma);
    printf("sigma          = ");
    for (int i = 0; i <= BCH_T; i++) printf("%02x ", sigma[i]);
    printf("\n");

    uint64_t fix = bch_syn_correction(error_codeword, sigma);
    printf("fix = %#llx\n", (unsigned long long)fix);
    printf("expected = %#llx\n", (unsigned long long)codeword);

    

    uint32_t found_og_data = bch_original_data(fix);
    printf("og_data = %#x\n", found_og_data);
    printf("expected = %#x\n", DATA);

    uint32_t decode = bch_decode(error_codeword);
    printf("direct decoded funtion to get data from errored codeword = %#x\n", decode);

    return 0;
}