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

constexpr int DGS_DISC_GAUSS_MAX_TABLE_SIZE_BYTES = (1 << 16);
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

long dgs_bern_mp_call(dgs_bern_mp_t *self) {
	self->prg->random_mpfr(self->tmp, 100); // Arbitrarally chose 100 bits as precision.

	if (mpfr_cmp(self->tmp, self->p)<0) {
		return 1;
	} else {
		return 0;
	}
}

void dgs_bern_mp_clear(dgs_bern_mp_t *self) {
	mpfr_clear(self->tmp);
	mpfr_clear(self->p);
	free(self);
}

static inline void _dgs_disc_gauss_mp_init_f(mpfr_t f, const mpfr_t sigma) {
	// we are assuming, mpfr_init was called on f already
	mpfr_set(f, sigma, MPFR_RNDN);
	mpfr_sqr(f, f, MPFR_RNDN);        // f = σ²
	mpfr_mul_ui(f, f, 2, MPFR_RNDN);  // f = 2 σ²
	mpfr_ui_div(f, 1, f, MPFR_RNDN);  // f = 1/(2 σ²)
	mpfr_neg(f, f, MPFR_RNDN);        // f = -1/(2 σ²)
}

static inline void _dgs_disc_gauss_mp_init_upper_bound(mpz_t upper_bound, mpz_t upper_bound_minus_one,
                                                       mpz_t two_upper_bound_minus_one, mpfr_t sigma, size_t tailcut) {
	mpfr_t tmp;
	mpfr_init2(tmp, mpfr_get_prec(sigma));
	mpz_init(upper_bound);
	mpz_init(upper_bound_minus_one);
	mpz_init(two_upper_bound_minus_one);
	mpfr_mul_ui(tmp, sigma, tailcut, MPFR_RNDN);        // tmp = σ·τ
	mpfr_add_ui(tmp, tmp, 1, MPFR_RNDN);                // tmp = σ·τ + 1
	mpfr_get_z(upper_bound, tmp, MPFR_RNDU);            // upper_bound = ⌈σ·τ + 1⌉
	mpz_sub_ui(upper_bound_minus_one, upper_bound, 1);  // upper_bound - 1 = ⌈σ·τ⌉
	mpz_mul_ui(two_upper_bound_minus_one, upper_bound, 2);
	mpz_sub_ui(two_upper_bound_minus_one, two_upper_bound_minus_one, 1);  // 2·upper_bound - 1
	mpfr_clear(tmp);
}

static inline void _dgs_disc_gauss_mp_init_rho(dgs_disc_gauss_mp_t *self, const mpfr_prec_t prec) {
	self->rho = (mpfr_t *)malloc(sizeof(mpfr_t) * mpz_get_ui(self->two_upper_bound_minus_one));
	if (!self->rho) {
		dgs_disc_gauss_mp_clear(self);
		throw std::bad_alloc();
	}

	mpfr_t x_;
	mpfr_init2(x_, prec);
	long absmax = mpz_get_ui(self->upper_bound) - 1;
	for (long x = -absmax; x <= absmax; x++) {
		mpfr_set_si(x_, x, MPFR_RNDN);
		mpfr_sub(x_, x_, self->c_r, MPFR_RNDN);
		mpfr_sqr(x_, x_, MPFR_RNDN);
		mpfr_mul(x_, x_, self->f, MPFR_RNDN);
		mpfr_exp(x_, x_, MPFR_RNDN);
		mpfr_init2(self->rho[x + absmax], prec);
		mpfr_set(self->rho[x + absmax], x_, MPFR_RNDN);
	}
	mpfr_clear(x_);
}


static inline long _dgs_disc_gauss_mp_min_in_rho(dgs_disc_gauss_mp_t *self, long range) {
	long mi   = 0;
	mpfr_t *m = &(self->rho[mi]);
	for (long x = 1; x < range; ++x) {
		if (mpfr_cmp(self->rho[x], *m) < 0) {
			mi = x;
			m  = &(self->rho[mi]);
		}
	}
	return mi;
}

static inline long _dgs_disc_gauss_mp_max_in_rho(dgs_disc_gauss_mp_t *self, long range) {
	long mi   = 0;
	mpfr_t *m = &(self->rho[mi]);
	for (long x = 1; x < range; ++x) {
		if (mpfr_cmp(self->rho[x], *m) > 0) {
			mi = x;
			m  = &(self->rho[mi]);
		}
	}
	return mi;
}


void dgs_disc_gauss_mp_call_alias(mpz_t rop, dgs_disc_gauss_mp_t *self) {
	//mpz_urandomm(rop, state, self->two_upper_bound_minus_one);
	self->prg->random_mpz(rop, self->two_upper_bound_minus_one);
	unsigned long x = mpz_get_ui(rop);
	if (self->bias[x]) {
		if (!dgs_bern_mp_call(self->bias[x])) {
			if (!self->alias[x]) {
				free(self);
				throw "bias initialized but no alias!";
			}
			mpz_set(rop, self->alias[x]);
		}
	}
	mpz_sub(rop, rop, self->upper_bound_minus_one);
	mpz_add(rop, rop, self->c_z);
}

void dgs_disc_gauss_mp_call_convolution(mpz_t rop, dgs_disc_gauss_mp_t *self) {
	mpz_set_si(rop, 0);
	for (int i = 0; i < self->n_coefficients; ++i) {
		self->base_sampler->call(self->x, self->base_sampler);
		mpz_mul(self->x, self->x, self->coefficients[i]);
		mpz_add(rop, rop, self->x);
	}

	// adjust center if necessary
	if (self->coset_sampler) {
		self->coset_sampler->call(self->x, self->coset_sampler);
		mpz_add(rop, rop, self->x);
	}

	mpz_add(rop, rop, self->c_z);
}



dgs_disc_gauss_mp_t *dgs_disc_gauss_mp_init(const mpfr_t sigma, const mpfr_t c, size_t tau,
                                            emp::PRG *prg, dgs_disc_gauss_alg_t algorithm) {
	if (mpfr_cmp_ui(sigma,0)<= 0)
		throw std::invalid_argument("Sigma must be > 0");
	if (tau == 0)
		throw std::invalid_argument("tau must be > 0");

	mpfr_prec_t prec                  = mpfr_get_prec(sigma);
	if (mpfr_get_prec(c) > prec) prec = mpfr_get_prec(c);

	dgs_disc_gauss_mp_t *self = (dgs_disc_gauss_mp_t *)calloc(sizeof(dgs_disc_gauss_mp_t), 1);
	if (!self) {
		throw std::bad_alloc();
	}

	mpz_init(self->x);
	mpz_init(self->x2);
	mpz_init(self->k);
	mpfr_init2(self->y, prec);
	mpfr_init2(self->z, prec);
	mpfr_init2(self->f, prec);

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

	switch (algorithm) {


	case DGS_DISC_GAUSS_ALIAS: {
		_dgs_disc_gauss_mp_init_upper_bound(self->upper_bound, self->upper_bound_minus_one, self->two_upper_bound_minus_one,
		                                    self->sigma, self->tau);
		_dgs_disc_gauss_mp_init_f(self->f, sigma);

		self->call = dgs_disc_gauss_mp_call_alias;
		if (mpz_cmp_ui(self->two_upper_bound_minus_one, ULONG_MAX / sizeof(mpfr_t)) > 0) {
			dgs_disc_gauss_mp_clear(self);
			throw std::overflow_error("integer overflow");
		}
		// we'll use the big table
		_dgs_disc_gauss_mp_init_rho(self, prec);

		// convert rho to probabilities
		mpfr_set_d(self->y, 0.0, MPFR_RNDN);
		mpfr_set_d(self->z, 1.0, MPFR_RNDN);
		long range = mpz_get_ui(self->two_upper_bound_minus_one);
		for (long x = 0; x < range; x++) { mpfr_add(self->y, self->y, self->rho[x], MPFR_RNDN); }
		mpfr_div(self->y, self->z, self->y, MPFR_RNDN);

		for (long x = 0; x < range; x++) { mpfr_mul(self->rho[x], self->rho[x], self->y, MPFR_RNDN); }

		// compute bias and alias
		self->alias = (mpz_t *)calloc(range, sizeof(mpz_t));
		if (!self->alias) {
			dgs_disc_gauss_mp_clear(self);
			throw std::bad_alloc();
		}

		self->bias = (dgs_bern_mp_t **)calloc(range, sizeof(dgs_bern_mp_t *));
		if (!self->bias) {
			dgs_disc_gauss_mp_clear(self);
			throw std::bad_alloc();
		}

		//~ // simple robin hood strategy approximates good alias
		//~ // this precomputation takes ~n^2, but could be reduced by
		//~ // using better data structures to compute min and max
		//~ // (instead of just linear search each time)
		mpfr_set_d(self->y, (double)range, MPFR_RNDN);
		mpfr_div(self->y, self->z, self->y, MPFR_RNDD);  // self->y = avg

		long low = _dgs_disc_gauss_mp_min_in_rho(self, range);
		long high;
		mpfr_sub(self->z, self->y, self->rho[low], MPFR_RNDD);  // z = avg - rho[low]

		mpfr_t p;
		mpfr_init2(p, prec);

		// we must be done after range rounds
		int n = 0;
		while (mpfr_cmp_d(self->z, 0) > 0 && n < range) {
			high = _dgs_disc_gauss_mp_max_in_rho(self, range);
			mpfr_mul_z(p, self->rho[low], self->two_upper_bound_minus_one, MPFR_RNDN);

			self->bias[low] = dgs_bern_mp_init(p, self->prg);
			mpz_init(self->alias[low]);
			mpz_set_ui(self->alias[low], high);
			mpfr_sub(self->rho[high], self->rho[high], self->z, MPFR_RNDU);
			mpfr_set(self->rho[low], self->y, MPFR_RNDU);

			low = _dgs_disc_gauss_mp_min_in_rho(self, range);
			mpfr_sub(self->z, self->y, self->rho[low], MPFR_RNDD);  // z = avg - rho[low]
			++n;
		}

		mpfr_clear(p);
		break;
	}

	case DGS_DISC_GAUSS_CONVOLUTION: {
		self->call = dgs_disc_gauss_mp_call_convolution;

		//~ double sigma1 = sigma;
		//~ double coset_sigma = sqrt(2)*eta;
		mpfr_t eta, sigma1, coset_sigma;
		mpfr_inits2(prec, eta, sigma1, coset_sigma, NULL);

		//~ double eta = sqrt((p+1)*log(2))/pi;
		mpfr_const_log2(eta, MPFR_RNDN);
		mpfr_const_pi(self->y, MPFR_RNDN);
		mpfr_mul_si(eta, eta, prec + 1, MPFR_RNDN);
		mpfr_sqrt(eta, eta, MPFR_RNDN);
		mpfr_div(eta, eta, self->y, MPFR_RNDN);

		mpfr_set(sigma1, sigma, MPFR_RNDN);
		mpfr_set_si(coset_sigma, 2, MPFR_RNDN);
		mpfr_sqrt(coset_sigma, coset_sigma, MPFR_RNDN);
		mpfr_mul(coset_sigma, coset_sigma, eta, MPFR_RNDN);

		if (!mpfr_zero_p(self->c_r)) {
			// we might need to adjust the center
			//~ sigma1 = sqrt(sigma*sigma - coset_sigma*coset_sigma);
			mpfr_mul(sigma1, sigma1, sigma1, MPFR_RNDN);

			mpfr_set(self->z, coset_sigma, MPFR_RNDN);
			mpfr_mul(self->z, self->z, self->z, MPFR_RNDN);

			mpfr_sub(sigma1, sigma1, self->z, MPFR_RNDN);
			mpfr_sqrt(sigma1, sigma1, MPFR_RNDN);
		}

		long table_size     = 2 * ceil(mpfr_get_ui(sigma1, MPFR_RNDU) * tau) * (sizeof(dgs_bern_mp_t) + sizeof(mpz_t));
		int recursion_level = 0;
		// for computing the recursion level, we can probably get away with doubles:
		double current_sigma = mpfr_get_d(sigma1, MPFR_RNDN);
		long z1              = 1;
		long z2              = 1;

		// compute recursion level for convolution
		while (table_size > (DGS_DISC_GAUSS_MAX_TABLE_SIZE_BYTES >> 1)) {
			recursion_level++;
			z1 = floor(sqrt(current_sigma / (mpfr_get_d(eta, MPFR_RNDN) * 2)));
			if (z1 == 0) {
				dgs_disc_gauss_mp_clear(self);
				throw "MAX_TABLE_SIZE too small to store alias sampler!";
			}
			z2 = (z1 > 1) ? z1 - 1 : 1;
			current_sigma /= (sqrt(z1 * z1 + z2 * z2));
			table_size = 2 * ceil(current_sigma * tau) * (sizeof(dgs_bern_mp_t) + sizeof(mpz_t));
		}

		self->n_coefficients = 1 << recursion_level;
		self->coefficients   = (mpz_t *)calloc(self->n_coefficients, sizeof(mpz_t));
		if (!self->coefficients) {
			dgs_disc_gauss_mp_clear(self);
			throw std::bad_alloc();
		}
		for (int i = 0; i < self->n_coefficients; ++i) {
			//~ self->coefficients[i] = 1;
			mpz_init(self->coefficients[i]);
			mpz_set_si(self->coefficients[i], 1);
		}

		// if there is no convolution, we simply forward to alias and
		// so we won't need adjustment of sigma, even if sampling off center
		//~ current_sigma = (recursion_level == 0)? sigma : sigma1;

		if (recursion_level == 0) {
			mpfr_set(self->y, sigma, MPFR_RNDN);
		} else {
			mpfr_set(self->y, sigma1, MPFR_RNDN);
		}

		// redo above computation to store coefficients
		for (int i = 0; i < recursion_level; ++i) {
			//~ z1 = floor(sqrt(current_sigma/(eta*2)));
			mpfr_set(self->z, self->y, MPFR_RNDN);
			mpfr_div(self->z, self->z, eta, MPFR_RNDN);
			mpfr_div_si(self->z, self->z, 2, MPFR_RNDN);
			mpfr_sqrt(self->z, self->z, MPFR_RNDN);
			mpfr_get_z(self->x, self->z, MPFR_RNDZ);

			//~ z2 = (z1 > 1) ? z1 - 1 : 1;
			int tmp = (mpz_cmp_si(self->x, 1) > 0);
			mpz_sub_ui(self->x2, self->x, tmp);

			// we unroll the recursion on the coefficients on the fly
			// so we don't have to use recursion during the call
			int off = (1 << recursion_level - i - 1);
			for (int j = 0; j < (1 << i); ++j) {
				for (int k = 0; k < off; ++k) {
					//~ self->coefficients[2*j*off + k] *= z1;
					mpz_mul(self->coefficients[2 * j * off + k], self->coefficients[2 * j * off + k], self->x);
				}
			}

			for (int j = 0; j < (1 << i); ++j) {
				for (int k = 0; k < off; ++k) {
					//~ self->coefficients[(2*j + 1)*off + k] *= z2;
					mpz_mul(self->coefficients[(2 * j + 1) * off + k], self->coefficients[(2 * j + 1) * off + k], self->x2);
				}
			}

			//~ current_sigma /= (sqrt(z1*z1 + z2*z2));
			mpz_mul(self->x, self->x, self->x);
			mpz_mul(self->x2, self->x2, self->x2);
			mpfr_set_z(self->z, self->x, MPFR_RNDN);
			mpfr_add_z(self->z, self->z, self->x2, MPFR_RNDN);
			mpfr_sqrt(self->z, self->z, MPFR_RNDN);
			mpfr_div(self->y, self->y, self->z, MPFR_RNDN);
		}

		//~ double base_c = self->c_r;
		mpfr_set(self->z, self->c_r, MPFR_RNDN);
		self->coset_sampler = NULL;
		//~ if (recursion_level > 0 && fabs(self->c_r) > 0) {
		if (recursion_level > 0 && !mpfr_zero_p(self->c_r)) {
			// we'll need to adjust the center
			//~ base_c = 0.0;
			mpfr_set_zero(self->z, 0);
			self->coset_sampler = dgs_disc_gauss_mp_init(coset_sigma, self->c_r, tau, self->prg, DGS_DISC_GAUSS_ALIAS);
		}

		self->base_sampler = dgs_disc_gauss_mp_init(self->y, self->z, tau, self->prg, DGS_DISC_GAUSS_ALIAS);

		mpfr_clears(eta, sigma1, coset_sigma, (mpfr_ptr)NULL);

		break;
	}

	default: free(self); throw std::invalid_argument("unknown algorithm");//, algorithm);
	}
	return self;
}

dgs_disc_gauss_mp_t *dgs_disc_gauss_mp_init(const mpfr_t sigma, const mpfr_t c, size_t tau, emp::PRG *prg) {
	return dgs_disc_gauss_mp_init(sigma, c, tau, prg, DGS_DISC_GAUSS_CONVOLUTION);
                                        
}


void dgs_disc_gauss_mp_clear(dgs_disc_gauss_mp_t *self) {
	mpfr_clear(self->sigma);
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
		unsigned long range                                           = mpz_get_ui(self->two_upper_bound_minus_one);

		for (unsigned long x = 0; x < range; x++) { mpfr_clear(self->rho[x]); }
		free(self->rho);
	}

	if (self->alias) {
		for (unsigned long x = 0; x < mpz_get_ui(self->two_upper_bound_minus_one); x++) {
			if (self->alias[x]) mpz_clear(self->alias[x]);
		}
		free(self->alias);
	}

	if (self->bias) {
		for (unsigned long x = 0; x < mpz_get_ui(self->two_upper_bound_minus_one); x++) {
			if (self->bias[x]) dgs_bern_mp_clear(self->bias[x]);
		}
		free(self->bias);
	}

	if (self->upper_bound) mpz_clear(self->upper_bound);
	if (self->upper_bound_minus_one) mpz_clear(self->upper_bound_minus_one);
	if (self->two_upper_bound_minus_one) mpz_clear(self->two_upper_bound_minus_one);

	if (self->base_sampler) { dgs_disc_gauss_mp_clear(self->base_sampler); }
	if (self->coefficients) {
		for (int i = 0; i < self->n_coefficients; ++i) {
			if (self->coefficients[i]) mpz_clear(self->coefficients[i]);
		}
		free(self->coefficients);
	}
	if (self->coset_sampler) { dgs_disc_gauss_mp_clear(self->coset_sampler); }

	free(self);
}



#endif // DGS_H__
