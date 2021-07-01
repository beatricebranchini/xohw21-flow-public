/******************************************
*MIT License
*
# *Copyright (c) [2020] [Beatrice Branchini, Luisa Cicolini, Giulia Gerometta]
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

extern "C" {
	#include "gap_affine/affine_wavefront_align.h"
}
#include <cstdlib>
#include <string>
#include <vector>
#include <iostream>
#include <chrono>
#include <string>
#include <string.h>
#include <fstream>
#include <ap_int.h>

#define TEXT_SIZE 64
#define PATTERN_SIZE 64
#define NUM 100
#define UNROLL_MAX 128
#define PORT_WIDTH 512
#define WAVEFRONT_SIZE (TEXT_SIZE*4)

#define AFFINE_WAVEFRONT_OFFSET_NULL (-10)

#define NOW std::chrono::high_resolution_clock::now();

void wfa(ap_uint<512>* text_input, int* textLength_input, ap_uint<512>* pattern_input, int* patternLength_input,
		int match, int mismatch, int gap_opening, int gap_extension, int* scores, int num_couples);
//		, int M[][WAVEFRONT_SIZE], int I[][WAVEFRONT_SIZE], int D[][WAVEFRONT_SIZE]);

int main(int argc, char const *argv[]) {

	int num_couples = NUM;

	// Pattern & Text
	std::vector<std::string> input_pattern(num_couples);
	std::vector<std::string> input_text(num_couples);

	// int seed = time(NULL);
	int seed=1624374933;

	std::cout<<"Seed: "<<seed<<std::endl;

	char alphabet[4] = {'A', 'C', 'G', 'T'};

	srand(seed);

	char pattern[PATTERN_SIZE*NUM];
	char text[TEXT_SIZE*NUM];

	int iter=0;

	for (int i=0; i<num_couples; i++) {
		for (int j=0; j<TEXT_SIZE; j++) {
			input_text[i].push_back(alphabet[rand()%4]);
			text[iter]=input_text[i][j];
			iter++;
		}
	}

	iter=0;
	for (int i=0; i<num_couples; i++) {
		for (int j=0; j<PATTERN_SIZE; j++) {
			input_pattern[i].push_back(alphabet[rand()%4]);
			pattern[iter]=input_pattern[i][j];
			iter++;
		}
	}

	int patternLength[NUM];
	int textLength[NUM];

	for (int i=0; i<num_couples; i++) {
		patternLength[i] = PATTERN_SIZE;
		textLength[i] = TEXT_SIZE;
	}

    // Set penalties
	affine_penalties_t affine_penalties = {
		.match = 0,
		.mismatch = 3,
		.gap_opening = 5,
		.gap_extension = 1,
	};

	int scores[NUM];

	std::vector<edit_cigar_t> result_cigars_check(input_text.size());

	int smarco_scores[NUM];

	FILE *f_smarco;

	f_smarco = fopen ("smarco.txt", "w");

	int M[WAVEFRONT_SIZE][WAVEFRONT_SIZE];// = (int*)malloc(sizeof(int)*TEXT_SIZE*4);
	int I[WAVEFRONT_SIZE][WAVEFRONT_SIZE];// = (int*)malloc(sizeof(int)*TEXT_SIZE*4);
	int D[WAVEFRONT_SIZE][WAVEFRONT_SIZE];// = (int*)malloc(sizeof(int)*TEXT_SIZE*4);

	for(int i = 0; i < WAVEFRONT_SIZE; i++){
		for(int j = 0; j < WAVEFRONT_SIZE; j++){
			M[i][j]=-10;
			I[i][j]=-10;
			D[i][j]=-10;
		}
	}

	ap_uint<PORT_WIDTH> *p_pattern, *p_text;
	p_pattern=(ap_uint<PORT_WIDTH>*)pattern;
	p_text=(ap_uint<PORT_WIDTH>*)text;

	wfa(p_text, textLength, p_pattern, patternLength, affine_penalties.match, affine_penalties.mismatch, affine_penalties.gap_opening, affine_penalties.gap_extension, scores, num_couples);//, M, I, D);

	std::chrono::high_resolution_clock::time_point start = NOW;

	for(unsigned i = 0 ; i < num_couples; i++){
		// Allocate MM
		mm_allocator_t* const mm_allocator = mm_allocator_new(BUFFER_SIZE_8M);
		// Init Affine-WFA
		affine_wavefronts_t* affine_wavefronts = affine_wavefronts_new_complete(
			input_pattern[i].length(),input_text[i].length(),&affine_penalties,NULL,mm_allocator);
		// Align
		affine_wavefronts_align(affine_wavefronts,input_pattern[i].c_str(),input_pattern[i].length(),input_text[i].c_str(),input_text[i].length(), &smarco_scores[i]);

		int idx_lo = affine_wavefronts->mwavefronts[smarco_scores[i]]->lo;
		int idx_hi = affine_wavefronts->mwavefronts[smarco_scores[i]]->hi;

		for (int k=idx_lo; k<=idx_hi; k++){
			if (affine_wavefronts->mwavefronts[smarco_scores[i]]->offsets[k] != AFFINE_WAVEFRONT_OFFSET_NULL)
				fprintf(f_smarco, "%3d", affine_wavefronts->mwavefronts[smarco_scores[i]]->offsets[k]);
		}

		fprintf(f_smarco, "\n");

		result_cigars_check[i] = affine_wavefronts->edit_cigar;
		result_cigars_check[i].operations = (char*)malloc(result_cigars_check[i].end_offset-result_cigars_check[i].begin_offset);
		strncpy(result_cigars_check[i].operations, affine_wavefronts->edit_cigar.operations+result_cigars_check[i].begin_offset,result_cigars_check[i].end_offset-result_cigars_check[i].begin_offset);

		//score_smarco = edit_cigar_score_edit(&affine_wavefronts->edit_cigar);

		affine_wavefronts_delete(affine_wavefronts);
		mm_allocator_delete(mm_allocator);

	}

	std::chrono::high_resolution_clock::time_point end = NOW;
	std::chrono::duration<double> time = std::chrono::duration_cast<std::chrono::duration<double>>(end-start);

	fclose(f_smarco);

	std::ifstream f_fpga("fpga_matrices.txt");
	std::ifstream f_sw("smarco.txt");

	std::string tmp_smarco;
	std::string tmp_fpga;

	int correct_wf = 0;
	int wrong_wf = 0;

	for(int i=0; i<num_couples; i++){
		std::getline(f_sw, tmp_smarco);
		std::getline(f_fpga, tmp_fpga);

		if (tmp_smarco.compare(tmp_fpga) == 0)
		    correct_wf++;
		else
			wrong_wf++;

		tmp_smarco.clear();
		tmp_fpga.clear();
	}

	bool test_score = true;

	for (int i=0; i<num_couples; i++){
		if (scores[i]!=smarco_scores[i]) test_score=false;
	}

	std::cout<<"Correct wfs: "<<correct_wf<<std::endl;
	std::cout<<"Wrong wfs: "<<wrong_wf<<std::endl;
	if (test_score) std::cout<<"All scores correct"<<std::endl;
	else std::cout<<"Test failed"<<std::endl;

	std::cout<< "Alignment took: " << time.count() << " seconds" <<std::endl;

	return 0;
}
