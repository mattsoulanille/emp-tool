#include "emp-tool/emp-tool.h"
#include <iostream>
#include "emp-tool/utils/dgs.h"
using namespace std;
using namespace emp;

int main() {
	PRG prg;//using a random seed

	mpfr_t sigma, c;
	
	mpfr_init(sigma);
	mpfr_init(c);
	mpfr_set_d(sigma, 4.0, MPFR_RNDN);
	mpfr_set_d(c, 0.0, MPFR_RNDN);

	size_t tau = 8;
	dgs_disc_gauss_mp_t *dgs = dgs_disc_gauss_mp_init(sigma, c, tau, &prg);

	mpz_t rop;
	mpz_init(rop);
	cout << "Samples: ";
	for (int i = 0; i < 10; i++) {
		dgs->call(rop, dgs);
		cout << mpz_get_si(rop) << ", ";
	}
	cout << endl;

	return 0;
	
}
