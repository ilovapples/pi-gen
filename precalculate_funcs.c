#include "gmp/gmp.h"
#include <stdlib.h>
#include <pthread.h>
#include <stdio.h>

typedef struct {
	unsigned long n_iterations;
	mpz_t *ret;
} precalc_thread_args_t;

void* precalc_p(void* args) {
	precalc_thread_args_t* data = (precalc_thread_args_t*) args;
	unsigned long n_iterations =  data->n_iterations;

	mpz_t *pval_out = (mpz_t*) malloc(sizeof(mpz_t) * n_iterations);
	
	if (pval_out == NULL) {
		printf("Failed to allocate memory for precalculating P(1,n).\nTerminating...\n");
		exit(1);
	}

	printf("Calculating P(1,n) for n in [1,%lu]...\033[E", n_iterations);
	fflush(stdout);
	for (unsigned long i=0U; i < n_iterations; ++i) {
		// mpz_init_set_si(pval_out[i], -(6*i-1));
		mpz_init_set_si(pval_out[i], -1);
		if (i == 0) {
			mpz_set_si(pval_out[i], 1);
			continue;
		}
		mpz_t tmp0,tmp1,tmp2;
		mpz_init_set_ui(tmp0, 6);
		mpz_init_set_ui(tmp1, i);
		mpz_addmul(pval_out[i], tmp0, tmp1); // should do -1+6*i (which = 6*i-1)
		mpz_neg(pval_out[i], pval_out[i]);

		mpz_init_set_si(tmp2, -5);
		mpz_addmul(tmp2, tmp0, tmp1); // should do -5+6*i (which = 6*i-5)
		mpz_mul(pval_out[i], pval_out[i], tmp2);

		mpz_set_ui(tmp0, 2);
		mpz_init_set_si(tmp2, -1);
		mpz_addmul(tmp2, tmp0, tmp1); // should do -1+2*i
		mpz_mul(pval_out[i], pval_out[i], tmp2);
		
		// mpz_init_set_si(tmp1, (2*i-1));
		// mpz_init_set_si(tmp2, (6*i-5));
		// mpz_mul(pval_out[i], pval_out[i], tmp1);
		// mpz_mul(pval_out[i], pval_out[i], tmp2);

		mpz_mul(pval_out[i], pval_out[i], pval_out[i-1]);

		mpz_clear(tmp0); mpz_clear(tmp1); mpz_clear(tmp2);
	}
	// printf("Finished calculating P_VALS.\n");
	printf("\033[1;41fDONE\033[E");

	data->ret = pval_out;
	
	pthread_exit(NULL);
}
void* precalc_q(void* args) {
	precalc_thread_args_t* data = (precalc_thread_args_t*) args;
	unsigned long n_iterations = data->n_iterations;

	mpz_t *qval_out = (mpz_t*) malloc(sizeof(mpz_t) * n_iterations);

	printf("Calculating Q(1,n) for n in [1,%lu]...\033[E", n_iterations);
	fflush(stdout);
	for (unsigned long i=0U; i < n_iterations; ++i) {
		mpz_init(qval_out[i]);
		if (i == 0) {
			mpz_set_si(qval_out[i], 1);
			continue;
		}
		mpz_t some_big_number;
		mpz_init_set_str(some_big_number, "10939058860032000", 10);
		mpz_t j_cubed; mpz_init_set_ui(j_cubed, i);
		mpz_pow_ui(j_cubed, j_cubed, 3);
		mpz_mul(qval_out[i], some_big_number, j_cubed);

		mpz_mul(qval_out[i], qval_out[i], qval_out[i-1]);

		mpz_clear(some_big_number); mpz_clear(j_cubed);
	}
	// printf("Finished calculating Q_VALS.\n");
	printf("\033[2;41fDONE\033[E");

	data->ret = qval_out;

	pthread_exit(NULL);
}

// binary splitting techniques from https://en.wikipedia.org/wiki/Chudnovsky_algorithm
void precalc_vals(mpq_t* sval_out, unsigned long n_iterations) {
    
    // precalculation of P_VALS (the values are offset one down, since P_VALS[0] is P(1,1))
    // mpz_t *P_VALS = (mpz_t*) malloc(sizeof(mpz_t) * n_iterations);
    // if (P_VALS == NULL) { // If memory cannot be allocated for P_VALS (displays error message and exits):
    //     printf("Failed to allocate memory for precalculating P(1,n) vals.");
    //     exit(1);
    // } else { // If memory CAN be allocated for P_VALS (precalculates all the P_VALS):

	const unsigned long APPR_CORRECT_DIGITS = (unsigned long) (n_iterations * 14.18165);

	precalc_thread_args_t p_args = { .n_iterations = n_iterations, .ret = NULL };
	precalc_thread_args_t q_args = { .n_iterations = n_iterations, .ret = NULL };
	
	pthread_t P_thread, Q_thread;
	pthread_create(&P_thread, NULL, precalc_p, (void*) &p_args);
	pthread_create(&Q_thread, NULL, precalc_q, (void*) &q_args);

	pthread_join(P_thread, NULL);
	pthread_join(Q_thread, NULL);

	mpz_t *P_VALS = p_args.ret;
	mpz_t *Q_VALS = q_args.ret;
	
	// gmp_printf("sixth element of P_VALS: %Zd\n", P_VALS[5]);
	// gmp_printf("sixth element of Q_VALS: %Zd\n", Q_VALS[5]);
		
    // precalculation of S_VALS (the values are offset one down, since S_VALS[0] is S(1,1) = 0)
    if (sval_out == NULL) { // If memory cannot be allocated for S_VALS (displays error message and exits)
		printf("Failed to allocate memory (outside the scope of this function) for precalculating S(1,n)");
		exit(1);
    }
    
    printf("Calculating S(1,n) for n in [1,%lu]...\033[E", n_iterations);
    fflush(stdout);
    for (unsigned long i=0L; i < n_iterations; ++i) {
		mpq_init(sval_out[i]);
        mpq_set_ui(sval_out[i], 0, 1);
        if (i == 0) {
            continue;
        }
        // calculate S(1,i)
        mpz_t numerator; mpz_init_set(numerator, P_VALS[i]);
        mpz_t tmp0; mpz_init_set_ui(tmp0, 13591409);
        mpz_t tmp1; mpz_init_set_ui(tmp1, 545140134);
        mpz_t tmp2; mpz_init_set_ui(tmp2, i);
        mpz_addmul(tmp0, tmp1, tmp2);
        mpz_mul(numerator, numerator, tmp0);
        
        mpq_set_num(sval_out[i], numerator);
        mpq_set_den(sval_out[i], Q_VALS[i]);
        mpq_canonicalize(sval_out[i]);
        
        if (i > 1) {
            mpq_add(sval_out[i], sval_out[i], sval_out[i-1]);
        }

        mpz_clear(tmp0); mpz_clear(tmp1); mpz_clear(tmp2); 
        mpz_clear(numerator); 
        if (i%200 == 0) {
            printf("iteration complete - S(1,%lu)\n\
\0337\033[;f\033[2KShould calculate around %lu correct digits.\n\033[2K\n\0338", 
				i+1, APPR_CORRECT_DIGITS);
        }
    }
    printf("Finished calculating S(1,n) for n in [1,%lu]\n", n_iterations);

	free(P_VALS);
	free(Q_VALS);
}
