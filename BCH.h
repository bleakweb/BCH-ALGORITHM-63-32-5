
#ifndef BCH_H
#define BCH_H

#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>

// BCH(63, 32) Parameters 
#define BCH_N 63      
#define BCH_K 32      
#define BCH_T 5       
#define BCH_PARITY  27
#define GF_PRIM_POLY  0x43



void bch_init(void);

uint64_t bch_encode(uint32_t data);

void bch_compute_syndromes(uint64_t codeword, uint8_t syndromes[2 * BCH_T]);

void bch_syn_error_locator(const uint8_t syndromes[2 * BCH_T], uint8_t sigma[BCH_T + 1]);

uint64_t bch_syn_correction(uint64_t codeword, const uint8_t sigma[BCH_T + 1]);

uint32_t bch_original_data(uint64_t corrected_codeword);

uint32_t bch_decode(uint64_t codeword);





#endif

