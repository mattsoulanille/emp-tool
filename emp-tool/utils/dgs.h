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

#ifndef DGS_H__
#define DGS_H__
#include "emp-tool/utils/constants.h"
#include "emp-tool/utils/prg.h"
#include <gmp.h>
#include <math.h>
#include <mpfr.h>
#include <stdexcept>
#include "emp-tool/utils/dgs_structs.h"

dgs_bern_mp_t* dgs_bern_mp_init(mpfr_t p, emp::PRG *prg) {
	/* we allow 0 and 1 here for low precision */
	assert((mpfr_cmp_d(p, 0.0) >= 0) && (mpfr_cmp_d(p, 1.0) <= 0));
	
	dgs_bern_mp_t *self = (dgs_bern_mp_t*) malloc(sizeof(dgs_bern_mp_t));
	if (!self) throw std::bad_alloc();
	self->prg = prg;

	mpfr_init2(self->p, mpfr_get_prec(p));
	mpfr_set(self->p, p, MPFR_RNDN);
	mpfr_init2(self->tmp, mpfr_get_prec(p));
	return self;
}

static inline void _dgs_disc_gauss_mp_init_upper_bound(mpz_t upper_bound,
                                                       mpz_t upper_bound_minus_one,
                                                       mpz_t two_upper_bound_minus_one,
                                                       mpfr_t sigma, size_t tailcut) {
	mpfr_t tmp;
	mpfr_init2(tmp, mpfr_get_prec(sigma));
	mpz_init(upper_bound);
	mpz_init(upper_bound_minus_one);
	mpz_init(two_upper_bound_minus_one);
	mpfr_mul_ui(tmp, sigma, tailcut, MPFR_RNDN); // tmp = σ·τ
	mpfr_add_ui(tmp, tmp, 1, MPFR_RNDN); // tmp = σ·τ + 1
	mpfr_get_z(upper_bound, tmp, MPFR_RNDU); // upper_bound = ⌈σ·τ + 1⌉
	mpz_sub_ui(upper_bound_minus_one, upper_bound, 1); // upper_bound - 1 = ⌈σ·τ⌉
	mpz_mul_ui(two_upper_bound_minus_one, upper_bound, 2);
	mpz_sub_ui(two_upper_bound_minus_one, two_upper_bound_minus_one, 1); // 2·upper_bound - 1
	mpfr_clear(tmp);
}

const int DGS_BERN_EXP_ALLOC_BLOCK_SIZE = 16;

dgs_bern_exp_mp_t* dgs_bern_exp_mp_init(mpfr_t f, size_t l, emp::PRG *prg) {
	dgs_bern_exp_mp_t *self = (dgs_bern_exp_mp_t *)malloc(sizeof(dgs_bern_exp_mp_t));
	if (!self) throw std::bad_alloc();

	self->prg = prg;
	/* l == 0, means we use the precision of f to decide l */
	if (l == 0)
		l = SIZE_MAX;

	self->l = DGS_BERN_EXP_ALLOC_BLOCK_SIZE;
	self->p = (mpfr_t*)malloc(sizeof(mpfr_t)*self->l);
	if (!self->p) throw std::bad_alloc();
	self->B = (dgs_bern_mp_t**)malloc(sizeof(dgs_bern_mp_t)*self->l);
	if (!self->B) throw std::bad_alloc();

	mpfr_t tmp, tmp2;
	mpfr_init2(tmp2, mpfr_get_prec(f));
	mpfr_init2(tmp, mpfr_get_prec(f));
	mpfr_set(tmp, f, MPFR_RNDN); // f
	mpfr_pow_si(tmp, tmp, -1, MPFR_RNDN); // 1/f
	mpfr_neg(tmp, tmp, MPFR_RNDN); // -1/f

	for(size_t i=0; i<l; i++) {
		mpfr_exp(tmp2, tmp, MPFR_RNDN);
		if (mpfr_zero_p(tmp2)) {
			self->l = i;
			break;
		}
		if (i%DGS_BERN_EXP_ALLOC_BLOCK_SIZE == 0 && i!=0) {
			self->l += DGS_BERN_EXP_ALLOC_BLOCK_SIZE;
			self->l = (l>self->l) ? self->l : l;
			self->p = (mpfr_t*) realloc(self->p, sizeof(mpfr_t)*self->l);
			if(!self->p) throw std::bad_alloc();
			self->B = (dgs_bern_mp_t **) realloc(self->B, sizeof(dgs_bern_exp_mp_t)*self->l);
			if(!self->B) throw std::bad_alloc();
		}

		mpfr_init2(self->p[i], mpfr_get_prec(f));
		mpfr_set(self->p[i], tmp2, MPFR_RNDN);
		self->B[i] = dgs_bern_mp_init(self->p[i], self->prg);

		mpfr_mul_ui(tmp, tmp, 2, MPFR_RNDN);
	}
	if (l < self->l)
		self->l = l;
	mpfr_clear(tmp);
	mpfr_clear(tmp2);
	return self;
}

static inline void _dgs_disc_gauss_mp_init_bexp(dgs_disc_gauss_mp_t *self, mpfr_t sigma, mpz_t upper_bound) {
	mpfr_init2(self->f, mpfr_get_prec(sigma));
	mpfr_set(self->f, sigma, MPFR_RNDN); // f = σ
	mpfr_sqr(self->f, self->f, MPFR_RNDN); // f = σ²
	mpfr_mul_ui(self->f, self->f, 2, MPFR_RNDN); // f = 2 σ²
	size_t l = 2*mpz_sizeinbase(upper_bound, 2);
	self->Bexp = dgs_bern_exp_mp_init(self->f, l, self->prg);
}

long dgs_bern_mp_call(dgs_bern_mp_t *self) {
	self->prg->random_mpfr(self->tmp, 100); // Arbitrarally chose 100 bits as precision.

	if (mpfr_cmp(self->tmp, self->p)<0) {
		return 1;
	} else {
		return 0;
	}
}

long dgs_bern_exp_mp_call(dgs_bern_exp_mp_t *self, mpz_t x) {//, gmp_randstate_t state) {
	assert(mpz_sgn(x) >= 0);
	long int start = (mpz_sizeinbase(x, 2) < self->l) ? mpz_sizeinbase(x, 2) : self->l;

	for(long int i=start-1; i>=0; i--) {
		if (mpz_tstbit(x, i)) {
			if (dgs_bern_mp_call(self->B[i]) == 0) {
				return 0;
			}
		}
	}
	return 1;
}

void dgs_disc_gauss_mp_call_uniform_logtable(mpz_t rop, dgs_disc_gauss_mp_t *self) {//, gmp_randstate_t state) {
	do {
		self->prg->random_mpz(self->x, self->two_upper_bound_minus_one);
		mpz_sub(self->x, self->x, self->upper_bound_minus_one);
		mpz_mul(self->x2, self->x, self->x);
	} while (dgs_bern_exp_mp_call(self->Bexp, self->x2) == 0);
	mpz_set(rop, self->x);
	mpz_add(rop, rop, self->c_z);
}

void dgs_bern_uniform_clear(dgs_bern_uniform_t *self) {
	mpz_clear(self->tmp);
	free(self);
}

void dgs_bern_mp_clear(dgs_bern_mp_t *self) {
	mpfr_clear(self->tmp);
	mpfr_clear(self->p);
	free(self);
}

void dgs_bern_exp_mp_clear(dgs_bern_exp_mp_t *self) {
	if(!self)
		return;

	for(size_t i=0; i<self->l; i++) {
		mpfr_clear(self->p[i]);
		dgs_bern_mp_clear(self->B[i]);
	}
	if(self->p)
		free(self->p);
	if(self->B)
		free(self->B);
	free(self);
}

void dgs_disc_gauss_mp_clear(dgs_disc_gauss_mp_t *self) {
	mpfr_clear(self->sigma);
	if (self->B) dgs_bern_uniform_clear(self->B);
	if (self->Bexp) dgs_bern_exp_mp_clear(self->Bexp);
	mpz_clear(self->x);
	mpz_clear(self->x2);
	mpz_clear(self->k);
	mpfr_clear(self->y);
	mpfr_clear(self->f);
	mpfr_clear(self->z);
	mpfr_clear(self->c);
	mpfr_clear(self->c_r);
	mpz_clear(self->y_z);
	mpz_clear(self->c_z);
	if (self->rho) {
		unsigned long range = mpz_get_ui(self->two_upper_bound_minus_one);
		for(unsigned long x=0; x<range; x++) {
			mpfr_clear(self->rho[x]);
		}
		free(self->rho);
	}
  
	if (self->alias) {
		for(unsigned long x=0; x<mpz_get_ui(self->two_upper_bound_minus_one); x++) {
			if (self->alias[x])
				mpz_clear(self->alias[x]);
		}
		free(self->alias);
	}
  
	if (self->upper_bound)
		mpz_clear(self->upper_bound);
	if (self->upper_bound_minus_one)
		mpz_clear(self->upper_bound_minus_one);
	if (self->two_upper_bound_minus_one)
		mpz_clear(self->two_upper_bound_minus_one);

	free(self);
}



dgs_disc_gauss_mp_t *dgs_disc_gauss_mp_init(const mpfr_t sigma, const mpfr_t c, size_t tau, emp::PRG *prg) {
	if (mpfr_cmp_ui(sigma,0)<= 0)
		throw std::invalid_argument("Sigma must be > 0");
	if (tau == 0)
		throw std::invalid_argument("tau must be > 0");

	mpfr_prec_t prec = mpfr_get_prec(sigma);
	if (mpfr_get_prec(c) > prec)
		prec = mpfr_get_prec(c);

	dgs_disc_gauss_mp_t *self = (dgs_disc_gauss_mp_t*) calloc(sizeof(dgs_disc_gauss_mp_t),1);
	if (!self) {
		throw std::bad_alloc();
	}

	mpz_init(self->x);
	mpz_init(self->x2);
	mpz_init(self->k);
	mpfr_init2(self->y, prec);
	mpfr_init2(self->z, prec);

	mpfr_init2(self->sigma, prec);
	mpfr_set(self->sigma, sigma, MPFR_RNDN);

	mpfr_init2(self->c, prec);
	mpfr_set(self->c, c, MPFR_RNDN);
	mpz_init(self->c_z);
	mpfr_get_z(self->c_z, c, MPFR_RNDN);
	mpfr_init2(self->c_r, prec);
	mpfr_sub_z(self->c_r, self->c, self->c_z, MPFR_RNDN);

	self->tau = tau;

	self->prg = prg;

	//self->algorithm = algorithm; // old code <<

	self->call = dgs_disc_gauss_mp_call_uniform_logtable;
	_dgs_disc_gauss_mp_init_upper_bound(self->upper_bound,
	                                    self->upper_bound_minus_one,
	                                    self->two_upper_bound_minus_one,
	                                    self->sigma, self->tau);

	if (!mpfr_zero_p(self->c_r)) {
		dgs_disc_gauss_mp_clear(self);
		throw std::invalid_argument("algorithm DGS_DISC_GAUSS_UNIFORM_LOGTABLE requires c%1 == 0");
	}

	_dgs_disc_gauss_mp_init_bexp(self, self->sigma, self->upper_bound);

	return self;
}





#endif // DGS_H__
