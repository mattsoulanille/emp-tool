/******************************************************************************                           
 *                                                                                                         
 *                        DGS - Discrete Gaussian Samplers                                                 
 *                                                                                                         
 * Copyright (c) 2014, Martin Albrecht  <martinralbrecht+dgs@googlemail.com>                               
 * All rights reserved.                                                                                    
 *                                                                                                         
 * Redistribution and use in source and binary forms, with or without                                      
 * modification, are permitted provided that the following conditions are met:                             
 *                                                                                                         
 * 1. Redistributions of source code must retain the above copyright notice, this                          
 *    list of conditions and the following disclaimer.                                                     
 * 2. Redistributions in binary form must reproduce the above copyright notice,                            
 *    this list of conditions and the following disclaimer in the documentation                            
 *    and/or other materials provided with the distribution.                                               
 *                                                                                                         
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"                             
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE                               
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE                          
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE                             
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL                              
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR                              
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER                              
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,                           
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE                           
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                                    
 *                                                                                                         
 * The views and conclusions contained in the software and documentation are                               
 * those of the authors and should not be interpreted as representing official                             
 * policies, either expressed or implied, of the FreeBSD Project.                                          
 ******************************************************************************/

#ifndef DGS_STRUCTS_H__
#define DGS_STRUCTS_H__
#include "emp-tool/utils/prg.h"
#include <gmp.h>
#include <math.h>
#include <mpfr.h>



typedef struct {

	/**
	   Number of bits we sample in each go.
	*/

	size_t   length;

	/**
	   Number of bits we have consumed yet.
	*/

	size_t   count;

	/**
	   We sample to this ``mpz_t``.
	*/

	mpz_t    tmp;

	/**
	   We store the pool of random bits here.
	*/

	unsigned long pool;
} dgs_bern_uniform_t;

typedef struct {
	/**                                                                                                 
        Return 1 with probability `p`                                                                    
	*/
	mpfr_t p;
	mpfr_t tmp; //< used internally
	emp::PRG *prg; // pseudorandom number generator
} dgs_bern_mp_t;

typedef struct {
	/**
	   We support positive `x` up to `2^l-1`
	*/

	size_t l;

	/**
	   Probabilities for Bernoulli sub-samplers.
	*/

	mpfr_t *p;

	/**
	   Bernoulli sub-samplers.
	*/

	dgs_bern_mp_t **B;

	emp::PRG *prg;

} dgs_bern_exp_mp_t;

typedef struct _dgs_disc_gauss_mp_t {

	/**
	   The width paramter `σ`, i.e. samples are accepted with probability
	   proportional to `\exp(-(x-c)²/(2σ²))`
	*/

	mpfr_t sigma;

	/**
	   The mean of the distribution `c`. The value of `c` does not have to be an
	   integer. However, some algorithms only support integer-valued `c`.
	*/

	mpfr_t c;

	mpfr_t c_r; //< `c_r := c % 1`
	mpz_t c_z;  //< c_z := c - (c_r)

	/**
	   Cutoff `τ`, samples outside the range `(⌊c⌉-⌈στ⌉,...,⌊c⌉+⌈στ⌉)` are
	   considered to have probability zero. This bound applies to algorithms
	   which sample from the uniform distribution.
	*/

	size_t tau;

	//dgs_disc_gauss_alg_t algorithm; //<  which algorithm to use

	/**
	   The Pseudorandom Number Generator. Source of randomness.
	*/
	emp::PRG* prg;

	
	/**
	   We use a uniform Bernoulli to decide signs.
	*/

	dgs_bern_uniform_t *B;

	/**
	   To realise rejection sampling, we call `B_{exp(-(x·x)/(2σ²))}` and accept
	   if it returns 1.

	   Used when ``DGS_DISC_GAUSS_UNIFORM_LOGTABLE`` or
	   ``DGS_DISC_GAUSS_SIGMA2_LOGTABLE`` is set.
	*/

	dgs_bern_exp_mp_t *Bexp;


	/**
	   Return an ``mpz_t`` sampled from this sampler

	   :param rop: target value.
	   :param self: discrete Gaussian sampler.
	   :param state: entropy pool.

	*/

	void (*call)(mpz_t rop, struct _dgs_disc_gauss_mp_t *self);//, gmp_randstate_t state);

	/**
	   We sample ``x`` with ``abs(x) < upper_bound`` in
	   ``DGS_DISC_GAUSS_UNIFORM_ONLINE``, ``DGS_DISC_GAUSS_UNIFORM_TABLE`` and
	   ``DGS_DISC_GAUSS_UNIFORM_LOGTABLE``.
	*/

	mpz_t upper_bound;

	/**
	   We sample ``x`` with ``abs(x) <= upper_bound - 1`` in
	   ``DGS_DISC_GAUSS_UNIFORM_ONLINE``, ``DGS_DISC_GAUSS_UNIFORM_TABLE`` and
	   ``DGS_DISC_GAUSS_UNIFORM_LOGTABLE``.
	*/

	mpz_t upper_bound_minus_one;

	/**
	   There are ``2*upper_bound -1`` elements in the range
	   ``-upper_bound+1,...,upper_bound-1``.
	*/

	mpz_t two_upper_bound_minus_one;

	/**
	   The multiplier `k` when we sample from `D_{k·σ₂,c}` in
	   ``DGS_DISC_GAUSS_SIGMA2_LOGTABLE``.
	*/

	mpz_t k;

	/**
	   Precomputed `-1/(2σ²)`.
	*/

	mpfr_t f;

	mpz_t x; //< space for temporary integer
	mpz_t y_z; //< space for temporary integer
	mpz_t x2; // space for temporary integer
	mpfr_t y; // space for temporary rational number
	mpfr_t z; // space for temporary rational number

	/**
	   Precomputed values for `exp(-(x-c)²/(2σ²))` in
	   ``DGS_DISC_GAUSS_UNIFORM_TABLE``
	*/

	mpfr_t *rho;
  
	/**
	 * Tables required for alias sampling.
	 */
   
	mpz_t* alias;
	//dgs_bern_mp_t** bias;

} dgs_disc_gauss_mp_t;

#endif
