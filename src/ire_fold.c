#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <stdbool.h>

#include "ire_fold.h"
#include "utils.h"

/*========================= Function Declarations =========================*/
IreStructure *predict_ire(const char *sequence);
void irestruct_destroy(IreStructure *ire_structure);
static IreStructure *irestruct_init(const char *sequence);
static bool worth_testing(const char *sequence);
static inline uint worth_testing_h(const char *sequence, int left, int right);
static inline uint predict_ire_h(IreStructure *ire_structure, int left, int right);
static inline uint get_max(IreStructure *ire, int left, int right);
static inline uint score_structure(const char *sequence, const char *structure);
static inline uint pair_is_gu(const char *sequence, int left, int right);
static inline uint check_pair(const char first_nucleotide, const char second_nucleotide);
static inline uint max3(uint a, uint b, uint c);


/*=========================  Library Functions   =========================*/
IreStructure *
predict_ire(const char *sequence)
{
	IreStructure *ire_structure = irestruct_init(sequence);
	if(!worth_testing(sequence)) {
		ire_structure->quality = 0;
		return ire_structure;
	}

	uint weight = predict_ire_h(ire_structure, 12, 19);
	if(weight != 0) {
		for(int i=0; i<13; i++) {
			if(ire_structure->structure[i] == '(') weight--;
		}
	}
	ire_structure->quality = weight;

	return ire_structure;
}


void
irestruct_destroy(IreStructure *ire_structure)
{
	free(ire_structure->structure);
	free(ire_structure->__struct);
	free(ire_structure);
}


/*========================= IRE Prediction Algorithm =========================*/
static inline uint
predict_ire_h(IreStructure *ire, int left, int right)
{
	if(left < 0 || right >= 32) {
		uint score=score_structure(ire->sequence,(const char *)ire->__struct);
		if(score >= ire->__bestscore) {
			strcpy(ire->structure, ire->__struct);
			ire->__bestscore = score;
		}
		return score;
	}
	
	return get_max(ire, left, right);
}


static inline uint
get_max(IreStructure *ire, int left, int right)
{
	uint score = check_pair(ire->sequence[left], ire->sequence[right]);

	if(score == 1) {
		ire->__struct[left] = '('; ire->__struct[right] = ')';
		score = predict_ire_h(ire, left-(1+(left==8)),right+1);
		ire->__struct[left] = '.'; ire->__struct[right] = '.';

	} else if(score == 2) {
		uint rscore = predict_ire_h(ire, left, right+1);
		uint lscore = predict_ire_h(ire, left-(1+(left==8)), right);
		ire->__struct[left] = '('; ire->__struct[right] = ')';
		uint mscore = predict_ire_h(ire, left-(1+(left==8)),right+1);
		ire->__struct[left] = '.'; ire->__struct[right] = '.';

		score = max3(rscore, lscore, mscore);

	} else {
		uint rscore = predict_ire_h(ire, left, right+1);
		uint lscore = predict_ire_h(ire, left-(1+(left==8)), right);
		uint mscore = predict_ire_h(ire, left-(1+(left==8)), right+1);

		score = max3(rscore, lscore, mscore);
	}

	return score;
}


static inline uint
score_structure(const char *sequence, const char *structure)
{
	uint score = 5;
	uint num_upper_pairs = 0;
	uint num_lower_pairs = 0;
	uint num_upper_bulge = 0;
	uint mismatch = 0;
	uint num_gu = 0;

	/* Check upper stem is valid, return 0 if not */
	for(int i=8; i<13; i++) {
		if(structure[i] == '(') num_upper_pairs++;
	}
	if(num_upper_pairs < 4) {
		return 0;
	}
	
	/* Check lower stem */
	for(int i=0; i<8; i++) {
		if(structure[i] == '(') num_lower_pairs++;
	}
	if(num_lower_pairs < 5) {
		return 0;
	}
	score += (num_upper_pairs + num_lower_pairs);

	/* Check pair types in upper stem */
	for(int i=12, j=19; i>5; i--) {
		if(i == 7) { continue; }

		if(structure[i] == '(' && structure[j] == ')') {
			num_gu += pair_is_gu(sequence, i, j);
		} else if(structure[i] == '.' && structure[j] == '.') {
			mismatch++;
		} else if(20 <= j && j <= 23 && num_upper_bulge < 1 && structure[j] == '.') {
			num_upper_bulge++;
			i++;
		} else {
			return 0;
		}
		j++;
	}

	if(mismatch && num_upper_bulge) {
		return 0;
	}
	if((num_gu > 2) || (mismatch > 1) || (num_upper_bulge > 1)) {
		return 0;
	}
	score -= num_gu;
	score -= num_upper_bulge;
	score -= mismatch;

	return score;
}


/*================== Helper Functions ==================*/
IreStructure *
irestruct_init(const char *sequence)
{
	IreStructure *ire_struct = s_malloc(sizeof *ire_struct);
	ire_struct->sequence  = sequence;
	ire_struct->structure = s_calloc(strlen(sequence) + 1, sizeof(char));
	ire_struct->__struct  = s_calloc(strlen(sequence) + 1, sizeof(char));
	ire_struct->__bestscore = 0;
	for(int i=0; i<(int)strlen(sequence); i++) {
		ire_struct->structure[i] = '.';
		ire_struct->__struct[i]  = '.';
	}

	return ire_struct;
}


static bool
worth_testing(const char *sequence)
{
	// Check upper stem to determine if worth predicting entire structure
	int num_pairs = worth_testing_h(sequence, 12, 19);
	return num_pairs > 4;
}


static inline uint
worth_testing_h(const char *sequence, int left, int right)
{
	if(left < 8 || right > 24) {
		return 0;
	}

	uint pair_type = check_pair(sequence[left], sequence[right]);

	uint rscore = worth_testing_h(sequence, left, right+1);
	uint lscore = worth_testing_h(sequence, left-1, right);
	uint mscore = worth_testing_h(sequence, left-1, right+1);

	uint score = max3(rscore, lscore, mscore);

	return score+(pair_type > 0);
}


static inline uint
pair_is_gu(const char *sequence, int left, int right)
{
	if(check_pair(sequence[left], sequence[right]) == 2) {
		return 1U;
	}
	return 0U;
}

static inline uint
check_pair(const char first_nucleotide, const char second_nucleotide)
{
	/* Define a mapping between nucleotide pairs and return values */
    const char         *pairs = "GC CG AT TA AU UA TG GT UG GU";
    const uint return_values[] = {1, 1, 1, 1, 1, 1, 2, 2, 2, 2};
	// const uint return_values[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    
    /* Iterate over pairs and compare with the input nucleotides */
    for (int i = 0; i < 10; ++i) {
        if (pairs[i * 3] == first_nucleotide && pairs[i * 3 + 1] == second_nucleotide) {
            return return_values[i];
        }
    }
    
    return 0U; /* If the pair is not found, return 0 */
}

static inline uint
max3(uint a, uint b, uint c)
{
    uint m = a;
    if (m < b) m = b;
    if (m < c) m = c;
    return m;
}
