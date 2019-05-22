//#include <cutil.h> //cuda toolkit below 5.0 support for CUDA_SAFE_CALL()
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <curand.h>
#include <curand_kernel.h>
#include <cuda_runtime.h>
//#include <helper_cuda.h>
#include <memory>

using namespace std;

// DEFINES 
#define NUM_BLOCKS 5*16//20*16 //5*16 //dimGrid //Keep numblocks a multiple of the #MP's of the GPU (8800GT=14MP)

//The register usage varies with platform. 64-bit Linux and 32.bit Windows XP have been tested.
#ifdef __linux__ //uses 25 registers per thread (64-bit)
#define NUM_THREADS_PER_BLOCK 320 //Keep above 192 to eliminate global memory access overhead However, keep low to allow enough registers per thread
#define NUM_THREADS 17920
#endif

#ifdef _WIN64
#define NUM_THREADS_PER_BLOCK 256  //256 //dimBlock
#define NUM_THREADS NUM_BLOCKS*NUM_THREADS_PER_BLOCK
#else //uses 26 registers per thread
#define NUM_THREADS_PER_BLOCK 288 //Keep above 192 to eliminate global memory access overhead However, keep low to allow enough registers per thread
#define NUM_THREADS 16128   
#endif


#define NUM_LAYER 6

#define NUMSTEPS_GPU       6000
#define PI                 3.141592654f
#define RPI                0.318309886f
#define MAX_LAYERS         100
#define STR_LEN            200
#define NORMAL             0                    // 1: normal, 0: oblique
#define NUM_OF_DETECTOR    (NORMAL ? 5:5) //(NORMAL ? 5:10)		//(NORMAL ? 4:9)       // normal: 4 fibers, oblique: 9 fibers
//#define ANGLE              (NORMAL ? 0:45)      // normal: 0 degree, oblique: 45 degree  
#define ANGLE              (NORMAL ? 0:0)      // normal: 0 degree, oblique: 0 degree  by CY

#define NAOfSource         (NORMAL ? 0.37:0.37)//(NORMAL ? 0.4:0.37)//(NORMAL ? 0.4:0.26)//(NORMAL ? 0.4:0.12)  // normal: 0.4, oblique; 0.22
#define NAOfDetector       (NORMAL ? 0.12:0.12)//(NORMAL ? 0.4:0.26)//(NORMAL ? 0.4:0.12)  // normal: 0.4, oblique; 0.22
#define n_detector         1.457//1.457//1.457//1.457 -fiber//1.61 //YU-modified
#define n_source           1.457//1.457//1.457//1.61 //YU-modified
#define illumination_r     0.075//0.075		//radius //Wang-modified //skin:0.025  IJV:0.075
#define collect_r          0.02//0.02//0.025//0.02			//radius //Wang-modified //skin:0.025  IJV:0.02
#define NUMBER_PHOTONS     1000000000//1000000000//50000000//400000000 -skin
#define NUMBER_SIMULATION  1//42//31//54//36  //IJV:36 skin:4

//#define WEIGHT 0.0001f
#define WEIGHTI 429497u //0xFFFFFFFFu*WEIGHT
#define CHANCE 0.1f

#define detected_temp_size 5000 //number of photon should be detected
#define SDS_detected_temp_size 10

// TYPEDEFS
typedef struct __align__(16)
{
	float z_min;		// Layer z_min [cm]
	float z_max;		// Layer z_max [cm]
	float mutr;			// Reciprocal mu_total [cm]
	float mua;			// Absorption coefficient [1/cm]
	float g;			// Anisotropy factor [-]
	float n;			// Refractive index [-]
}LayerStruct;

typedef struct __align__(16)
{
	float x;		// Global x coordinate [cm]
	float y;		// Global y coordinate [cm]
	float z;		// Global z coordinate [cm]
	float dx;		// (Global, normalized) x-direction
	float dy;		// (Global, normalized) y-direction
	float dz;		// (Global, normalized) z-direction
	float weight;			// Photon weight
	int layer;				// Current layer
	bool first_scatter; // flag of first scatter
	float layer_pathlength[NUM_LAYER]; // record the pathlength in each layer for detected photon
	int scatter_event; // record howmany time the photon had been scattered
}PhotonStruct;

typedef struct
{
	unsigned long long number_of_photons;
	unsigned int n_layers;
	float start_weight;
	LayerStruct* layers;
}SimulationStruct;

typedef struct
{
	float radius[NUM_OF_DETECTOR+1];
	float NA[NUM_OF_DETECTOR+1];
	float position[NUM_OF_DETECTOR+1];
	float angle[NUM_OF_DETECTOR+1];
	
	//bool photon_detected[NUM_OF_DETECTOR + 1]; // whether this fiber had detected photon
	int detected_photon_counter; // record how many photon had been detected by this fiber in one iteration, should not exceed SDS_detected_temp_size
	float data[SDS_detected_temp_size]; // the photon weight detected by this probe
	float layer_pathlength[SDS_detected_temp_size][NUM_LAYER]; // record the pathlength in each layer for detected photon
	int scatter_event[SDS_detected_temp_size]; // record howmany time the photon had been scattered
	int detected_SDS_number[SDS_detected_temp_size]; // record which SDS detected the photon
}Fibers;

typedef struct
{
	Fibers* f;
	PhotonStruct* p;					// Pointer to structure array containing all the photon data
	unsigned int* thread_active;		// Pointer to the array containing the thread active status
	unsigned long long* num_terminated_photons;	//Pointer to a scalar keeping track of the number of terminated photons
	curandState*  state;
}MemStruct;

typedef struct
{
	float* all;
	float* prob;
	float* cumf;
}G_Array;


__device__ __constant__ unsigned long long num_photons_dc[1];
__device__ __constant__ unsigned int n_layers_dc[1];
__device__ __constant__ float start_weight_dc[1];
__device__ __constant__ LayerStruct layers_dc[MAX_LAYERS];