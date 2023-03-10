/* (C) Copyright 1993-2006, Fred Hutchinson Cancer Research Center */
/* Use, modification or distribution of these programs is subject to */
/* the terms of the non-commercial licensing agreement in license.h. */

/* gcode.h: Definitions for the translation of nucleotides to amino acids */
/* Modified by: Bill Alford */
/* Change log information is at the end of the file. */

#ifndef __GCODE_H__
#define __GCODE_H__

typedef struct {
	char	*name;
	char	*code;
	char	*inits;
	} GeneticCode, *GeneticCodePtr;

#define CODON_LEN	3	/* No. of nucleotides per codon */

#define NUMBER_OF_GENETIC_CODES 17
EXTERN GeneticCode gcodes[NUMBER_OF_GENETIC_CODES]
#ifdef INIT
= {
{"Standard",
"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
"0001000000000000000100000000000000010000000000000000000000000000"},

{"Vertebrate Mitochondrial",
"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG",
"0000000000000000000000000000000011110000000000000001000000000000"},

{"Yeast Mitochondrial",
"FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
"0000000000000000000000000000000000110000000000000000000000000000"},

{"Mold Mitochondrial and Mycoplasma",
"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
"0011000000000000000100000000000011110000000000000001000000000000"},

{"Invertebrate Mitochondrial",
"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG",
"0001000000000000000000000000000011110000000000000001000000000000"},

{"Ciliate Nuclear",
"FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
"0000000000000000000000000000000000010000000000000000000000000000"},

{"Echinoderm Mitochondrial",
"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
"0000000000000000000000000000000000010000000000000001000000000000"},

{"Euplotid Nuclear" ,
"FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
"0000000000000000000000000000000000010000000000000000000000000000"},

{"Bacterial and Plant Plastid" ,
"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
"0001000000000000000100000000000011110000000000000001000000000000"},

{"Alternative Yeast Nuclear" ,
"FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
"0000000000000000000100000000000000010000000000000000000000000000"},

{"Ascidian Mitochondrial" ,
"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG",
"0000000000000000000000000000000000010000000000000000000000000000"},

{"Flatworm Mitochondrial" ,
"FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
"0000000000000000000000000000000000010000000000000000000000000000"},

{"Blepharisma Macronuclear" ,
"FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
"0000000000000000000000000000000000010000000000000000000000000000"},

{"Chlorophycean Mitochondrial" ,
"FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
"0000000000000000000000000000000000010000000000000000000000000000"},

{"Trematode Mitochondrial" ,
"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
"0000000000000000000000000000000000010000000000000001000000000000"},

{"Scenedesmus obliquus mitochondrial" ,
"FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
"0000000000000000000000000000000000010000000000000000000000000000"},

{"Thraustochytrium mitochondrial code" ,
"FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
"0000000000000000000000000000000010010000000000000001000000000000"}
}
#endif
	;


void init_gcode PROTO((GeneticCodePtr, unsigned char [64], unsigned char [64]));
unsigned char codon2aa PROTO((unsigned char *, unsigned, unsigned, unsigned));
void aa2codon PROTO((unsigned char *));

#endif /* !__GCODE_H__ */




/* Change log information follows. */
/* 
   Changes since version 3.6:
   10/14/04 fixed gcodes table (missing commas!)
   Changes since version 3.4:
   8/ 8/01  Updated from ncbi/data/gc.prt
   Changes since version 3.1:
   1/20/97  Added aa2codon()
 */
/*  Old codes
"Protozoan Mitochondrial",
"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
"0011000000000000000100000000000000010000000000000000000000000000",

"Plant Mitochondrial",
"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRWIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
"0000000000000000000000000000000000010000000000000000000000000000",
*/
