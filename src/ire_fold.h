#ifndef KATSS_IRE_FOLD_H
#define KATSS_IRE_FOLD_H

typedef unsigned int uint;

typedef struct IreStructure {
	const char *sequence;    /* Sequence of the associated IRE */
	char       *structure;   /* Predicted structure of the IRE */
	uint        quality;     /* Quality score of the predicted structure */

	char       *__struct;    /* Stores structure that is being processed */
	uint        __bestscore; /* Current best score in the preediction */
} IreStructure;

/**
 * @brief Predict the structure and quality of an IRE.
 * 
 * @param sequence The IRE sequence
 * @return IreStructure* 
 */
IreStructure *predict_ire(const char *sequence);


/**
 * @brief Remove an IreStructure an all of its associated resources.
 * 
 * @param ire_structure Pointer to the IreStructure struct
 */
void irestruct_destroy(IreStructure *ire_structure);
#endif // KATSS_IRE_FOLD_H
