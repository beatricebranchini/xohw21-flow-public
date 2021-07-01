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

extern "C" {
	#include "gap_affine/affine_wavefront_align.h"
}

#include <fstream>
#include "xcl2.hpp"
#include <algorithm>
#include <iostream>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <string>
#include <climits>
#include <time.h>
#include <chrono>

#define TEXT_SIZE 1024
#define PATTERN_SIZE 1024
#define NUM 100
#define PORT_WIDTH 512
#define MAX_NUM NUM
#define WAVEFRONT_SIZE (TEXT_SIZE*4)
#define AFFINE_WAVEFRONT_OFFSET_NULL (-10)

#define NUM_KERNEL 6

#define NOW std::chrono::high_resolution_clock::now();

#define MAX_HBM_BANKCOUNT 32
#define BANK_NAME(n) n | XCL_MEM_TOPOLOGY
const int bank[MAX_HBM_BANKCOUNT] = {
    BANK_NAME(0),  BANK_NAME(1),  BANK_NAME(2),  BANK_NAME(3),  BANK_NAME(4),
    BANK_NAME(5),  BANK_NAME(6),  BANK_NAME(7),  BANK_NAME(8),  BANK_NAME(9),
    BANK_NAME(10), BANK_NAME(11), BANK_NAME(12), BANK_NAME(13), BANK_NAME(14),
    BANK_NAME(15), BANK_NAME(16), BANK_NAME(17), BANK_NAME(18), BANK_NAME(19),
    BANK_NAME(20), BANK_NAME(21), BANK_NAME(22), BANK_NAME(23), BANK_NAME(24),
    BANK_NAME(25), BANK_NAME(26), BANK_NAME(27), BANK_NAME(28), BANK_NAME(29),
    BANK_NAME(30), BANK_NAME(31)};

int main(int argc, char *argv[]){
    
    std::string binaryFile;
	std::string readsPath;
	
    cl::Context context;
 
    cl::CommandQueue commands;
    std::vector< char, aligned_allocator<char> > pattern(MAX_NUM*PATTERN_SIZE);
    std::vector< char, aligned_allocator<char> > text(MAX_NUM*TEXT_SIZE);
    std::vector< int, aligned_allocator<int> > scores(MAX_NUM);
    std::vector< int, aligned_allocator<int> > patternLength(MAX_NUM);
    std::vector< int, aligned_allocator<int> > textLength(MAX_NUM);

    int num = MAX_NUM; //Number of patterns set to the maximum by default
    int seed = time(NULL);

    std::vector<std::string> pattern_test(NUM);
    std::vector<std::string> text_test(NUM);
	
    if (argc == 3) { //Input provided by file 

        binaryFile = argv[1];
        readsPath = argv[2];
        
        std::ifstream st(readsPath);
	    std::string temp1;
	    std::string temp2;
	    
        int idx = 0;
        while (st >> temp1 >> temp2){ //Reading input from file 
		    std::string temp;
		    temp = temp1.substr(1,std::string::npos);
		    pattern_test[idx] = temp;
		    temp = temp2.substr(1,std::string::npos);
            text_test[idx] = temp;
            idx++;
	    }

        num = idx;

        for (int i=0; i<NUM; i++){ //Storing lengths
            patternLength[i] = pattern_test[i].size();
            textLength[i] = text_test[i].size();
            if (i>=num){
                patternLength[i] = 0;
                textLength[i] = 0;
            }
        }

        int iter_pattern = 0;
        int iter_text = 0;
        for (int i=0; i<num; i++){ //Reorganizing inputs for the core
            
            for (int j=0; j<pattern_test[i].size(); j++){
                pattern[iter_pattern] = pattern_test[i][j];
                iter_pattern++;
            }
            
            for (int k=0; k<text_test[i].size(); k++){
                text[iter_text] = text_test[i][k];
                iter_text++;
            }
        }
    
    } else { //Random generation of inputs
 
        binaryFile = argv[1];

        std::cout<<"Seed: "<<seed<<std::endl;
	    srand(seed);

	    char alphabet[4] = {'A', 'C', 'G', 'T'};

	    int iter = 0; 
	    for (int i = 0; i < num; i++) {
		    for (int j = 0; j < TEXT_SIZE; j++) {
			    text_test[i].push_back(alphabet[rand()%4]);
			    text[iter]=text_test[i][j];
			    iter++;
		    }
	    }

	    iter = 0;
	    for (int i = 0; i < num; i++) {
		    for (int j = 0; j < PATTERN_SIZE; j++) {
			    pattern_test[i].push_back(alphabet[rand()%4]);
			    pattern[iter] = pattern_test[i][j];
			    iter++;
		    }
	    }

	    for (int i = 0; i < num; i++) {
		    patternLength[i] = PATTERN_SIZE;
		    textLength[i] = TEXT_SIZE;
	    }
    }

    // Set penalties
	affine_penalties_t affine_penalties = {
		.match = 0,
		.mismatch = 3,
		.gap_opening = 5,
		.gap_extension = 1,
	};

	if ((affine_penalties.mismatch-affine_penalties.gap_extension) != (affine_penalties.gap_opening-affine_penalties.mismatch)) {
		std::cout<<"Error: unsupported penalties!"<<std::endl;
		std::cout<<"Must satisfy: (mismatch - gap extension) == (gap opening - mismatch)"<<std::endl;
		std::cout<<"Test failed"<<std::endl;
		exit(1);
	}

	printf("PENALTIES INITIALIZED: %d, %d, %d, %d \n", affine_penalties.match, affine_penalties.mismatch, affine_penalties.gap_opening, affine_penalties.gap_extension);

    std::string krnl_name = "wfa";
    std::vector<cl::Kernel> krnls(NUM_KERNEL);

    // The get_xil_devices will return vector of Xilinx Devices
    auto devices = xcl::get_xil_devices();

    // read_binary_file() command will find the OpenCL binary file created using the
    // V++ compiler load into OpenCL Binary and return pointer to file buffer.
    auto fileBuf = xcl::read_binary_file(binaryFile);

    cl::Program::Binaries bins{{fileBuf.data(), fileBuf.size()}};
    int valid_device = 0;

    cl_int err;

    for (unsigned int i = 0; i < devices.size(); i++) {
        auto device = devices[i];
            // Creating Context and Command Queue for selected Device
        OCL_CHECK(err, context = cl::Context(device, NULL, NULL, NULL, &err));
        OCL_CHECK(err, commands = cl::CommandQueue(context, device,
                            CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE | CL_QUEUE_PROFILING_ENABLE, &err));

        std::cout << "Trying to program device[" << i 
                  << "]: " << device.getInfo<CL_DEVICE_NAME>() << std::endl;
     
        cl::Program program(context, {device}, bins, NULL, &err);
        
        if (err != CL_SUCCESS) {
            std::cout << "Failed to program device[" << i
                        << "] with xclbin file!\n";                      
        } else {
            std::cout << "Device[" << i << "]: program successful!\n";
            
            // Creating Kernel object using Compute unit names
            for (int i = 0; i < NUM_KERNEL; i++) {
                std::string cu_id = std::to_string(i + 1);
                std::string krnl_name_full = krnl_name + ":{" + "wfa_" + cu_id + "}";

                printf("Creating a kernel [%s] for CU(%d)\n", krnl_name_full.c_str(), i + 1);

                //Here Kernel object is created by specifying kernel name along with compute unit.
                //For such case, this kernel object can only access the specific Compute unit
                OCL_CHECK(err, krnls[i] = cl::Kernel(program, krnl_name_full.c_str(), &err));
            }

            valid_device++;
            break; // we break because we found a valid device
        }
        std::cout<<"dwvgae"<<std::endl;
    }

	std::cout<<"Kernel created"<<std::endl;
    
    if (valid_device == 0) {
        std::cout << "Failed to program any device found, exit!\n";
        exit(EXIT_FAILURE);
    }

    // Create device buffers
    std::vector<cl_mem_ext_ptr_t> pattern_buffer_ext(NUM_KERNEL);
    std::vector<cl_mem_ext_ptr_t> text_buffer_ext(NUM_KERNEL);
    std::vector<cl_mem_ext_ptr_t> scores_buffer_ext(NUM_KERNEL);
    std::vector<cl_mem_ext_ptr_t> textLength_buffer_ext(NUM_KERNEL);
    std::vector<cl_mem_ext_ptr_t> patternLength_buffer_ext(NUM_KERNEL);

    std::vector<cl::Buffer> pattern_buffer(NUM_KERNEL);
    std::vector<cl::Buffer> text_buffer(NUM_KERNEL);
    std::vector<cl::Buffer> scores_buffer(NUM_KERNEL);
    std::vector<cl::Buffer> textLength_buffer(NUM_KERNEL);
    std::vector<cl::Buffer> patternLength_buffer(NUM_KERNEL);


    for(int i = 0; i < NUM_KERNEL; i++) {

        text_buffer_ext[i].obj = text.data();
        text_buffer_ext[i].param = 0;
        text_buffer_ext[i].flags = bank[i*5];

        textLength_buffer_ext[i].obj = textLength.data();
        textLength_buffer_ext[i].param = 0;
        textLength_buffer_ext[i].flags = bank[i*5+1];
        
        pattern_buffer_ext[i].obj = pattern.data();
        pattern_buffer_ext[i].param = 0;
        pattern_buffer_ext[i].flags = bank[i*5+2];

        patternLength_buffer_ext[i].obj = patternLength.data();
        patternLength_buffer_ext[i].param = 0;
        patternLength_buffer_ext[i].flags = bank[i*5+3];

        scores_buffer_ext[i].obj = scores.data();
        scores_buffer_ext[i].param = 0;
        scores_buffer_ext[i].flags = bank[i*5+4];

    }

    for (int i = 0; i < NUM_KERNEL; i++) {
    	OCL_CHECK(err, text_buffer[i] = cl::Buffer(context, CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX |
                                    CL_MEM_USE_HOST_PTR, sizeof(char)*MAX_NUM*TEXT_SIZE, &text_buffer_ext[i], &err));
        OCL_CHECK(err, textLength_buffer[i] = cl::Buffer(context, CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX |
                                    CL_MEM_USE_HOST_PTR, sizeof(int)*MAX_NUM, &textLength_buffer_ext[i], &err));
        OCL_CHECK(err, pattern_buffer[i] = cl::Buffer(context, CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX |
                                    CL_MEM_USE_HOST_PTR, sizeof(char)*MAX_NUM*PATTERN_SIZE, &pattern_buffer_ext[i], &err));
        OCL_CHECK(err, patternLength_buffer[i] = cl::Buffer(context, CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX |
                                    CL_MEM_USE_HOST_PTR, sizeof(int)*MAX_NUM, &patternLength_buffer_ext[i], &err));
        OCL_CHECK(err, scores_buffer[i] = cl::Buffer(context, CL_MEM_READ_WRITE | CL_MEM_EXT_PTR_XILINX |
                                    CL_MEM_USE_HOST_PTR, sizeof(int)*MAX_NUM, &scores_buffer_ext[i], &err));
	}

	commands.finish();

    // Write our data set into device buffers  
     for(int i = 0; i < NUM_KERNEL; i++)
        err = commands.enqueueMigrateMemObjects({text_buffer[i], textLength_buffer[i], pattern_buffer[i], patternLength_buffer[i]}, 0);

    if (err != CL_SUCCESS) {
            printf("Error: Failed to write to device memory!\n");
            printf("Test failed\n");
            exit(1);
    }

	commands.finish();
    
    // Set the arguments to our compute kernel
    for (int i = 0; i < NUM_KERNEL; i++) {
        OCL_CHECK(err, err = krnls[i].setArg(0, text_buffer[i]));
        OCL_CHECK(err, err = krnls[i].setArg(1, textLength_buffer[i]));
        OCL_CHECK(err, err = krnls[i].setArg(2, pattern_buffer[i]));
        OCL_CHECK(err, err = krnls[i].setArg(3, patternLength_buffer[i]));
		OCL_CHECK(err, err = krnls[i].setArg(4, affine_penalties.match));
		OCL_CHECK(err, err = krnls[i].setArg(5, affine_penalties.mismatch));
		OCL_CHECK(err, err = krnls[i].setArg(6, affine_penalties.gap_opening));
		OCL_CHECK(err, err = krnls[i].setArg(7, affine_penalties.gap_extension));
        OCL_CHECK(err, err = krnls[i].setArg(8, scores_buffer[i]));
		OCL_CHECK(err, err = krnls[i].setArg(9, num));

        if (err != CL_SUCCESS) {
            printf("Error: Failed to set kernel arguments! %d\n", err);
            printf("Test failed\n");
            exit(1);
        }
    }

	commands.finish();

    std::chrono::high_resolution_clock::time_point start = NOW;
    
    // Execute the kernel over the entire range of our 1d input data set
    // using the maximum number of work group items for this device
    for (int i = 0; i < NUM_KERNEL; ++i)
        err |= commands.enqueueTask(krnls[i]);


    if (err) {
        printf("Error: Failed to execute kernel! %d\n", err);
        printf("Test failed\n");
        exit(1);
    }

    commands.finish();
    std::chrono::high_resolution_clock::time_point end = NOW;
	std::chrono::duration<double> time = std::chrono::duration_cast<std::chrono::duration<double>>(end-start)/NUM_KERNEL;
    
    // Read back the results from the device to verify the output
    for (int i = 0; i < NUM_KERNEL; ++i) {
        err = commands.enqueueMigrateMemObjects({scores_buffer[i]}, CL_MIGRATE_MEM_OBJECT_HOST);  
    }


    if (err != CL_SUCCESS) {
        printf("Error: Failed to read output array! %d\n", err);
        printf("Test failed\n");
        exit(1);
    }

	printf("HW time: %lf\n", time);

	//Checking the results 
    int sw_scores[MAX_NUM];

    start = NOW;

    for(unsigned i = 0 ; i < num; i++){
		// Allocate MM
		mm_allocator_t* const mm_allocator = mm_allocator_new(BUFFER_SIZE_8M);
		// Init Affine-WFA
		affine_wavefronts_t* affine_wavefronts = affine_wavefronts_new_complete(
			pattern_test[i].size(), text_test[i].length(), &affine_penalties, NULL, mm_allocator);
		// Align
		affine_wavefronts_align(affine_wavefronts, pattern_test[i].c_str(), pattern_test[i].length(),
			text_test[i].c_str(), text_test[i].length(), &sw_scores[i]);

		affine_wavefronts_delete(affine_wavefronts);
		mm_allocator_delete(mm_allocator);

	}

	end = NOW;
	time = std::chrono::duration_cast<std::chrono::duration<double>>(end-start);

	printf("SW time: %lf\n", time);

	bool test_score = true;

	for (int i=0; i<num; i++){
		if (scores[i]!=sw_scores[i]){
            printf("HW: %d, SW: %d\n", scores[i], sw_scores[i]);
            test_score=false;
        }
	}

	if (test_score) 
		std::cout<<"All scores correct"<<std::endl;
	else 
		std::cout<<"Test failed"<<std::endl;

}