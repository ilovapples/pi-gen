#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <time.h>

#include "precalculate_funcs.c"

// TODO: Implement rolling values to reduce space complexity from (3 * O(n)) to (3 * O(1)).
// 		 We're only using the last value of S_VALS, so the necessary P_VALS and S_VALS can be 
//       calculated in the calculation for S_VALS when they're needed, with 3 variables that
//       are overwritten with each new term so that we can still access the last term for the
//       overall product/summation.

void calc_final_value_to(mpf_t out, mpq_t *S_VALS, unsigned long n_iterations) {
	printf("Calculating pi with %lu iterations of a binary-split variation of the Chudnovsky algorithm.\n", n_iterations);

	printf("Putting together the terms to calculate pi...");

	// TERM 1 (entire numerator apparently)
	mpf_t term10, term11;
	mpf_init_set_ui(term10, 426880);
	mpf_init(term11); 
	mpf_sqrt_ui(term11, 10005);
	mpf_t numer; mpf_init(numer); 
	mpf_mul(numer, term10, term11);

	// TERM 2 (entire denominator i guess)
	mpq_t a_big_number; mpq_init(a_big_number);
	mpq_set_ui(a_big_number, 13591409, 1);
	mpq_t denom; mpq_init(denom);
	mpq_add(denom, a_big_number, S_VALS[n_iterations-1]);

	// put it all together
	mpf_t mpf_denom; mpf_init(mpf_denom);
	mpf_set_q(mpf_denom, denom);
	mpf_div(out, numer, mpf_denom);
	
	printf("DONE\n");
}

const int APPR_MAX_DIGITS = 200000;
const int IS_DEBUG = 1;

int main(int argc, char* argv[]) {
	const int APPR_SUFFICIENT_BITS = (int) (APPR_MAX_DIGITS * 3.5);

	mpf_set_default_prec(APPR_SUFFICIENT_BITS);
	
	mpf_t pi_get;
	mpf_init_set_ui(pi_get, 0);
	mpf_set_prec(pi_get, APPR_SUFFICIENT_BITS);

	unsigned long n_iterations = 30UL;

	if (argc > 1) {
		n_iterations = atoi(argv[1]);
	}

	printf("\033[;f\033[2J");

	clock_t start = clock();

	mpq_t *S_VALS = (mpq_t*) malloc(sizeof(mpq_t) * n_iterations);
	precalc_vals(S_VALS, n_iterations); // precalculate pvals and qvals for the algorithm

	calc_final_value_to(pi_get, S_VALS, n_iterations);	

	clock_t end = clock();
	double elapsed_secs = (double)(end - start) / CLOCKS_PER_SEC;

	mp_exp_t exp;
	char* significand = mpf_get_str(
		NULL,
		&exp,
		10,
		APPR_MAX_DIGITS,
		pi_get);
	
	// gmp_printf("significand: %.Ff", pi_get);
	printf("significand: \n%s\n", significand);

	
	FILE *fptr = fopen("check.txt", "r");
	unsigned long correct_to = 0;
	char check[1000000];
	fgets(check, sizeof(check), fptr);
	fclose(fptr);
	for (unsigned long i=0; i < strlen(significand); ++i) {
		if (significand[i] != check[i]) {
			correct_to = i - 1; // minus one to ignore the 3
			if (i >= strlen(check)) {
				printf("Calculated more digits of pi than the reference file has! You should consult a larger reference file to get the true accuracy.\n");
			} else {
				printf("This approximation of pi was accurate to %lu decimal places!\n", correct_to);
			}
			break;
		}
	}
	if (correct_to == 0) {
		printf("All %lu calculated digits are correct! Try a higher input to see how much higher you can go.", n_iterations);
	}
	
	fptr = fopen("out.txt", "w");
	fprintf(fptr, "%s", significand);
	fclose(fptr);

	printf("The entire calculation took %.9f seconds.", elapsed_secs);

	// free variables
	for (unsigned long i=0; i < n_iterations; ++i) {
		mpq_clear(S_VALS[i]);
	}
	free(S_VALS);

	return 0;
}
