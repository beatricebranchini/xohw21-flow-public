/******************************************
*MIT License
*
# *Copyright (c) [2020] [Beatrice Branchini, Luisa Cicolini, Giulia Gerometta, Marco Santambrogio]
*
*Permission is hereby granted, free of charge, to any person obtaining a copy
*of this software and associated documentation files (the "Software"), to deal
*in the Software without restriction, including without limitation the rights
*to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
*copies of the Software, and to permit persons to whom the Software is
*furnished to do so, subject to the following conditions:
*
*The above copyright notice and this permission notice shall be included in all
*copies or substantial portions of the Software.
*
*THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
*IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
*FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
*AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
*LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
*OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
*SOFTWARE.
******************************************/

#include <cstdlib>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <ap_int.h>

#define TEXT_SIZE 1024
#define PATTERN_SIZE 1024
#define NUM 100
#define WAVEFRONT_SIZE (TEXT_SIZE*4)
#define N_ELEM_BLOCK 64
#define CACHED_WAVEFRONTS 9

#define UNROLL_FACTOR 16

const unsigned int unroll_f = UNROLL_FACTOR;

//#define WAVEFRONT_DEPTH (TEXT_SIZE*2-1)
//#define WAVEFRONT_WIDTH (TEXT_SIZE*2+2)
#define UNROLL_MAX 128
#define PORT_WIDTH 512

#define AFFINE_WAVEFRONT_OFFSET_NULL (-10)

typedef struct {
  int match;              // (Penalty representation; usually M <= 0)
  int mismatch;           // (Penalty representation; usually X > 0)
  int gap_opening;        // (Penalty representation; usually O > 0)
  int gap_extension;      // (Penalty representation; usually E > 0)
} affine_penalties_t;

int compute_max(int n1, int n2, int n3, int n4) {
	int max;
	max = (n1 > n2) ? (n1 > n3 ? n1 : n3) : (n2 > n3 ? n2 : n3);
	if (max>n4) return max;
	else return n4;
}

int compute_max3(int n1, int n2, int n3) {
	return (n1 > n2) ? (n1 > n3 ? n1 : n3) : (n2 > n3 ? n2 : n3);
}

int compute_max2(int n1, int n2) {
	if (n1>n2) return n1;
	else return n2;
}

int compute_min(int n1, int n2, int n3, int n4) {
	int min;
	min = (n1 < n2) ? (n1 < n3 ? n1 : n3) : (n2 < n3 ? n2 : n3);
	if (min<n4) return min;
	else return n4;
}

template< typename t, unsigned int WAVE_SIZE >
void compute_wavefront_limit(t M[][WAVE_SIZE], int position, int* M_hi_sx, int* M_lo_sx){

	loop1_compute_wf_limit: for(int i = 0; i < WAVEFRONT_SIZE/2; i++){
#pragma HLS pipeline
		t wavefront_value = M[position%CACHED_WAVEFRONTS][WAVEFRONT_SIZE/2 + i];
		if( wavefront_value != AFFINE_WAVEFRONT_OFFSET_NULL)
			*M_hi_sx = i;
	}

	loop2_compute_wf_limit: for(int i = 0; i < WAVEFRONT_SIZE/2; i++){
#pragma HLS pipeline
		t wavefront_value = M[position%CACHED_WAVEFRONTS][WAVEFRONT_SIZE/2 - i];
		if( wavefront_value != AFFINE_WAVEFRONT_OFFSET_NULL)
			*M_lo_sx = -i;
	}
}

void align_m(short M[][WAVEFRONT_SIZE], int lo, int hi, int s, affine_penalties_t* affine_penalties){

	loop_align_m: for (int k = 0; k < WAVEFRONT_SIZE; k+=unroll_f) {

#pragma HLS pipeline
		for(int i = 0; i < unroll_f; i++){
#pragma HLS unroll
#pragma HLS dependence variable=M inter false
			if(WAVEFRONT_SIZE/2+lo<=k+i && k+i<=WAVEFRONT_SIZE/2+hi){
				M[s%CACHED_WAVEFRONTS][k + i] = 1 + M[(s-affine_penalties->mismatch)%CACHED_WAVEFRONTS][k+i];
			}
		}
	}

}

void align_dm(short M[][WAVEFRONT_SIZE], short D[][WAVEFRONT_SIZE], int lo, int hi, int s, affine_penalties_t* affine_penalties){

	loop_align_dm_0: for (int k = lo; k <= hi; k++) {
#pragma HLS pipeline
#pragma HLS dependence variable=M inter false
#pragma HLS dependence variable=D inter false
		D[s%CACHED_WAVEFRONTS][(WAVEFRONT_SIZE/2) + k] = compute_max2(
					M[(s-affine_penalties->gap_opening-affine_penalties->gap_extension)%CACHED_WAVEFRONTS][(WAVEFRONT_SIZE/2)+k+1],
					D[(s-affine_penalties->gap_extension)%CACHED_WAVEFRONTS][(WAVEFRONT_SIZE/2)+k+1]);
	}

	loop_align_dm_1: for (int k = lo; k <= hi; k++) {
#pragma HLS pipeline
#pragma HLS dependence variable=M inter false
#pragma HLS dependence variable=D inter false
		M[s%CACHED_WAVEFRONTS][(WAVEFRONT_SIZE/2) + k] = compute_max2(
					1 + M[(s-affine_penalties->mismatch)%CACHED_WAVEFRONTS][(WAVEFRONT_SIZE/2)+k],
					D[s%CACHED_WAVEFRONTS][(WAVEFRONT_SIZE/2)+k]);
	}
}

void align_im(short M[][WAVEFRONT_SIZE], short I[][WAVEFRONT_SIZE], int lo, int hi, int s, affine_penalties_t* affine_penalties){

	loop_align_im_0: for (int k = lo; k <= hi; k++) {
#pragma HLS pipeline
#pragma HLS dependence variable=M inter false
#pragma HLS dependence variable=I inter false
		I[s%CACHED_WAVEFRONTS][(WAVEFRONT_SIZE/2) + k] = 1 + compute_max2(
							M[(s-affine_penalties->gap_opening-affine_penalties->gap_extension)%CACHED_WAVEFRONTS][(WAVEFRONT_SIZE/2)+k-1],
							I[(s-affine_penalties->gap_extension)%CACHED_WAVEFRONTS][(WAVEFRONT_SIZE/2)+k-1]);
	}

	loop_align_im_1: for (int k = lo; k <= hi; k++) {
#pragma HLS pipeline
#pragma HLS dependence variable=M inter false
#pragma HLS dependence variable=I inter false
		M[s%CACHED_WAVEFRONTS][(WAVEFRONT_SIZE/2) + k] = compute_max2(
							1 + M[(s-affine_penalties->mismatch)%CACHED_WAVEFRONTS][(WAVEFRONT_SIZE/2)+k],
							I[s%CACHED_WAVEFRONTS][(WAVEFRONT_SIZE/2)+k]);
	}
}

void align_idm(short M[][WAVEFRONT_SIZE], short D[][WAVEFRONT_SIZE], short I[][WAVEFRONT_SIZE], int lo, int hi, int s, affine_penalties_t* affine_penalties){



	loop_align_idm_0: for (int k = 0; k < WAVEFRONT_SIZE; k+=unroll_f) {
#pragma HLS pipeline
		for(int i = 0; i < unroll_f; i++){
#pragma HLS unroll
#pragma HLS dependence variable=M inter false
#pragma HLS dependence variable=I inter false
			if(WAVEFRONT_SIZE/2+lo<=k+i && k+i<=WAVEFRONT_SIZE/2+hi){
				I[s%CACHED_WAVEFRONTS][k + i] = 1 + compute_max2(
										M[(s-affine_penalties->gap_opening-affine_penalties->gap_extension)%CACHED_WAVEFRONTS][i+k-1],
										I[(s-affine_penalties->gap_extension)%CACHED_WAVEFRONTS][i+k-1]);
			}
		}
	}

	loop_align_idm_1: for (int k = 0; k < WAVEFRONT_SIZE; k+=unroll_f) {
#pragma HLS pipeline
		for(int i = 0; i < unroll_f; i++){
#pragma HLS unroll
#pragma HLS dependence variable=M inter false
#pragma HLS dependence variable=D inter false
			if(WAVEFRONT_SIZE/2+lo<=k+i && k+i<=WAVEFRONT_SIZE/2+hi){
				D[s%CACHED_WAVEFRONTS][i + k] = compute_max2(
											M[(s-affine_penalties->gap_opening-affine_penalties->gap_extension)%CACHED_WAVEFRONTS][i+k+1],
											D[(s-affine_penalties->gap_extension)%CACHED_WAVEFRONTS][i+k+1]);
			}
		}
	}


	loop_align_idm_2: for (int k = 0; k < WAVEFRONT_SIZE; k+=unroll_f) {
#pragma HLS pipeline
		for(int i = 0; i < unroll_f; i++){
#pragma HLS dependence variable=M inter false
#pragma HLS dependence variable=I inter false
#pragma HLS dependence variable=D inter false
			if(WAVEFRONT_SIZE/2+lo<=k+i && k+i<=WAVEFRONT_SIZE/2+hi){
				M[s%CACHED_WAVEFRONTS][i + k] = compute_max3(
							1 + M[(s-affine_penalties->mismatch)%CACHED_WAVEFRONTS][i+k],
							I[s%CACHED_WAVEFRONTS][i+k],
							D[s%CACHED_WAVEFRONTS][i+k]);
			}
		}
	}
}

extern "C" {
void wfa(ap_uint<512>* text_input, int* textLength_input, ap_uint<512>* pattern_input, int* patternLength_input,
		int match, int mismatch, int gap_opening, int gap_extension, int* scores, int num_couples){

#pragma HLS INTERFACE m_axi port=text_input offset=slave bundle=gmem0 //depth=2
#pragma HLS INTERFACE m_axi port=textLength_input offset=slave bundle=gmem1 //depth=2
#pragma HLS INTERFACE m_axi port=pattern_input offset=slave bundle=gmem2 //depth=2
#pragma HLS INTERFACE m_axi port=patternLength_input offset=slave bundle=gmem3 //depth=2
#pragma HLS INTERFACE m_axi port=scores offset=slave bundle=gmem4 //depth=2

#pragma HLS INTERFACE s_axilite port=text_input bundle=control
#pragma HLS INTERFACE s_axilite port=textLength_input bundle=control
#pragma HLS INTERFACE s_axilite port=pattern_input bundle=control
#pragma HLS INTERFACE s_axilite port=patternLength_input bundle=control
#pragma HLS INTERFACE s_axilite port=match bundle=control
#pragma HLS INTERFACE s_axilite port=mismatch bundle=control
#pragma HLS INTERFACE s_axilite port=gap_opening bundle=control
#pragma HLS INTERFACE s_axilite port=gap_extension bundle=control
#pragma HLS INTERFACE s_axilite port=scores bundle=control
#pragma HLS INTERFACE s_axilite port=num_couples bundle=control

#pragma HLS INTERFACE s_axilite port=return bundle=control


	ap_uint<PORT_WIDTH> p_pattern[(TEXT_SIZE*NUM-1)/N_ELEM_BLOCK+1];
	ap_uint<PORT_WIDTH> p_text[(TEXT_SIZE*NUM-1)/N_ELEM_BLOCK+1];

	char *pattern = (char*)p_pattern;
	char *text = (char*)p_text;
	
	int textLength[NUM], patternLength[NUM];

	affine_penalties_t affine_penalties;

	int scores_local[NUM];

	int k, s, pos_s;
	int hi, lo, v, h;
	int M_lo_sx, M_hi_sx, M_lo_soe, M_hi_soe, I_lo_se, I_hi_se, D_lo_se, D_hi_se;

	int pattern_count = 0;
	int text_count = 0;

	short M[CACHED_WAVEFRONTS][WAVEFRONT_SIZE];
#pragma HLS ARRAY_PARTITION variable=M dim=1 complete
#pragma HLS ARRAY_PARTITION variable=M dim=2 factor=unroll_f cyclic

	short I[CACHED_WAVEFRONTS][WAVEFRONT_SIZE];
#pragma HLS ARRAY_PARTITION variable=I dim=1 complete
#pragma HLS ARRAY_PARTITION variable=I dim=2 factor=unroll_f cyclic

	short D[CACHED_WAVEFRONTS][WAVEFRONT_SIZE];
#pragma HLS ARRAY_PARTITION variable=D dim=1 complete
#pragma HLS ARRAY_PARTITION variable=D dim=2 factor=unroll_f cyclic

	//Initialization of wavefronts to the offset
	initialize_wavefronts:for(int i = 0; i < CACHED_WAVEFRONTS; i++){
		for(int j = 0; j < WAVEFRONT_SIZE; j++){
			M[i][j] = AFFINE_WAVEFRONT_OFFSET_NULL;
			I[i][j] = AFFINE_WAVEFRONT_OFFSET_NULL;
			D[i][j] = AFFINE_WAVEFRONT_OFFSET_NULL;
		}

	}

	//Local copies for the reads
	loop_text_pattern_copy: for (int i = 0; i < TEXT_SIZE*NUM/64; i++) {
#pragma HLS pipeline
		p_text[i] = text_input[i];
		p_pattern[i] = pattern_input[i];
	}

	//Local copies for lenghts of inputs 
	loop_data_copy: for (int i=0; i<NUM; i++) {
#pragma HLS pipeline
		textLength[i] = textLength_input[i];
		patternLength[i] = patternLength_input[i];
	}
	
	affine_penalties.match = match;
	affine_penalties.mismatch = mismatch;
	affine_penalties.gap_extension = gap_extension;
	affine_penalties.gap_opening = gap_opening;

	//Default initialization 
	M[0][WAVEFRONT_SIZE/2] = 0;
	s = 0; //Optimal alignment score, number of wavefront 

	//Loop iterating on the couples of reads
	loop_couples: for (unsigned n=0; n<num_couples; n++) {

		//Loop aligning the reads
		loop_alignment: while(true){

			//WF EXTEND - Algorithm 2 from the paper
			//Extension of the wavefront along consecutive matches
			wf_extend: for (k = (-patternLength[n]); k <= textLength[n]; k++){

				v = M[s%CACHED_WAVEFRONTS][WAVEFRONT_SIZE/2 + k] - k;
				h = M[s%CACHED_WAVEFRONTS][WAVEFRONT_SIZE/2 + k];

				if (v>=0 && h>=0) {

					//Support variables
					int idx_pattern = pattern_count + v;
					int idx_text = text_count + h;
					char char_pattern = pattern[idx_pattern];
					char char_text = text[idx_text];
					bool test_toast = pattern[idx_pattern] != text[idx_text];
					bool old_test_toast = 0;

					inner_wf_extend: while(v<patternLength[n] && h<textLength[n]){
#pragma HLS pipeline
						test_toast = pattern[idx_pattern] != text[idx_text];
						if(old_test_toast) {
							break;
						}
						idx_pattern++;
						idx_text++;
						v++;
						h++;
						old_test_toast = test_toast;
					}
					M[s%CACHED_WAVEFRONTS][WAVEFRONT_SIZE/2 + k] = h - old_test_toast;
				}
			}

			//If the alignment is complete, move to the next couple
			if (M[s%CACHED_WAVEFRONTS][WAVEFRONT_SIZE/2 + (textLength[n] - patternLength[n])] >= textLength[n]){

				scores_local[n]=s;
				text_count = text_count + textLength[n];
				pattern_count = pattern_count + patternLength[n];
 
 				for(int k = 0 ; k < CACHED_WAVEFRONTS; k++){
 					for (int j =  0; j < WAVEFRONT_SIZE; j++) {

 	 					M[k][j] = AFFINE_WAVEFRONT_OFFSET_NULL;
						I[k][j] = AFFINE_WAVEFRONT_OFFSET_NULL;
						D[k][j] = AFFINE_WAVEFRONT_OFFSET_NULL;

 					}
 				}

 				M[0%CACHED_WAVEFRONTS][WAVEFRONT_SIZE/2] = 0;
				s = 0;

				break;
			}

			s++;

			M[s%CACHED_WAVEFRONTS][WAVEFRONT_SIZE/2] = AFFINE_WAVEFRONT_OFFSET_NULL;
			scores_local[n] = s;

			// WF NEXT - Algorithm 3 from the paper
			// Computation of the next wavefront
			
			// Check null wavefronts
			int varM1=0;
			int varD=0;
			int varI=0;
			int varM2=0;

			// If any var == 0 the corresponding wavefront is null
			if ((s-affine_penalties.gap_opening-affine_penalties.gap_extension>=0) && 
				(M[(s-affine_penalties.gap_opening-affine_penalties.gap_extension)%CACHED_WAVEFRONTS][WAVEFRONT_SIZE/2]!=AFFINE_WAVEFRONT_OFFSET_NULL)){
				varM1=1;
			}
			if ((s-affine_penalties.gap_extension>=0) &&
				(D[(s-affine_penalties.gap_extension)%CACHED_WAVEFRONTS][WAVEFRONT_SIZE/2]!=AFFINE_WAVEFRONT_OFFSET_NULL)){
				varD=1;
			}
			if((s-affine_penalties.gap_extension>=0) &&
				(I[(s-affine_penalties.gap_extension)%CACHED_WAVEFRONTS][WAVEFRONT_SIZE/2]!=AFFINE_WAVEFRONT_OFFSET_NULL)){
				varI=1;
			}
			if((s-affine_penalties.mismatch>=0) &&
				(M[(s-affine_penalties.mismatch)%CACHED_WAVEFRONTS][WAVEFRONT_SIZE/2]!=AFFINE_WAVEFRONT_OFFSET_NULL)){
				varM2=1;
			}

			// Wavefront alignment procedure
			int tmp_mismatch = s-affine_penalties.mismatch;
			int tmp_gap_extension = s-affine_penalties.gap_extension;
			int tmp_gap_open = s-affine_penalties.gap_opening-affine_penalties.gap_extension;
			int control = varM1 + varD + varI + varM2;

			if ( ((tmp_mismatch >= 0) || (tmp_gap_extension >= 0) || (tmp_gap_open >= 0)) && ((control)>0)){

				// Default initialization 
				M_hi_sx = -1;
				M_lo_sx = 1;
				I_hi_se = -1;
				I_lo_se = 1;
				D_hi_se = -1;
				D_lo_se = 1;
				M_hi_soe = -1;
				M_lo_soe = 1;

				// Searching the upper and lower bounds of the wavefronts
				compute_wavefront_limit<short, WAVEFRONT_SIZE>(M, s-affine_penalties.mismatch, &M_hi_sx, &M_lo_sx);
				compute_wavefront_limit<short, WAVEFRONT_SIZE>(I, s-affine_penalties.gap_extension, &I_hi_se, &I_lo_se);
				compute_wavefront_limit<short, WAVEFRONT_SIZE>(D, s-affine_penalties.gap_extension, &D_hi_se, &D_lo_se);

				compute_wavefront_limit<short,WAVEFRONT_SIZE>(M, s-affine_penalties.gap_opening-affine_penalties.gap_extension, &M_hi_soe, &M_lo_soe);

				if (s-affine_penalties.mismatch < 0){
					M_hi_sx = -1;
					M_lo_sx = 1;
				}

				if (s-affine_penalties.gap_extension < 0){
					I_hi_se = -1;
					I_lo_se = 1;
					D_hi_se = -1;
					D_lo_se = 1;
				}

				if (s-affine_penalties.gap_opening-affine_penalties.gap_extension < 0){
					M_hi_soe = -1;
					M_lo_soe = 1;
				}
				

				hi = compute_max(M_hi_sx, M_hi_soe, I_hi_se, D_hi_se) + 1;
				lo = compute_min(M_lo_sx, M_lo_soe, I_lo_se, D_lo_se) - 1;

				if ((s-affine_penalties.gap_opening-affine_penalties.gap_extension>=0)&&(s-affine_penalties.gap_extension>=0)){
					varI=1;
					varD=1;
				} else {
					varI=0;
					varD=0;
				}

				int kernel;

				//Computation of values of the wavefronts 
				kernel = (varI << 1) | (varD);
				switch (kernel) {

					case 3: // 11b
						align_idm(M, D, I, lo, hi, s, &affine_penalties);
						break;

					case 2: // 10b
						align_im(M, I, lo, hi, s, &affine_penalties);
						break;

					case 1: // 01b
						align_dm(M, D, lo, hi, s, &affine_penalties);
						break;

					case 0: // 00b
						align_m(M, lo, hi, s, &affine_penalties);
						break;

				}
			}
		}

	}

	// Copying back the results
	result_copy: for (int i=0; i<NUM; i++) {
#pragma HLS pipeline
		scores[i]=scores_local[i];
	}
}
}