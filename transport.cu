#include "header.h"
//#include <helper_cuda.h>	//YU-modified
//#include <helper_string.h>  //YU-modified
//#include <helper_math.h>	//YU-modified
#include <float.h> //for FLT_MAX


int InitMemStructs(MemStruct* HostMem, MemStruct* DeviceMem, SimulationStruct* sim, char* fiber1_position); //Wang modified
void FreeMemStructs(MemStruct* HostMem, MemStruct* DeviceMem);
void FreeSimulationStruct(SimulationStruct* sim, int n_simulations);
__global__ void MCd(MemStruct DeviceMem, unsigned long long seed);
__global__ void LaunchPhoton_Global(MemStruct DeviceMem);
int InitDCMem(SimulationStruct* sim);
int Write_Simulation_Results(MemStruct* HostMem, SimulationStruct* sim, clock_t simulation_time);
int read_simulation_data(char* filename, SimulationStruct** simulations, int ignoreAdetection);
int interpret_arg(int argc, char* argv[], unsigned long long* seed, int* ignoreAdetection);

__global__ void MCd(MemStruct DeviceMem);
__device__ void LaunchPhoton(PhotonStruct* p, curandState *state);
__global__ void LaunchPhoton_Global(MemStruct DeviceMem, unsigned long long seed);
__device__ void Spin(PhotonStruct*, float, curandState* state);
__device__ unsigned int Reflect(PhotonStruct*, int, curandState* state);
__device__ unsigned int PhotonSurvive(PhotonStruct*, curandState* state);
__device__ void AtomicAddULL(unsigned long long* address, unsigned int add);
__device__ void detect(PhotonStruct* p, Fibers* f);
__device__ int binarySearch(float *data, float value);
void fiber_initialization(Fibers* f, float fiber1_position); //Wang modified
void output_fiber(SimulationStruct* sim, float* reflectance, char* output); //Wang modified
void calculate_reflectance(Fibers* f, float *result);
void input_g(int index, G_Array *g);
int InitG(G_Array* HostG, G_Array* DeviceG, int index);
void FreeG(G_Array* HostG, G_Array* DeviceG);

void output_A_rz(SimulationStruct* sim, unsigned long long *data, char* output);
void output_A0_z(SimulationStruct* sim, unsigned long long *data, char* output);


__device__ float rn_gen(curandState *s)
{
	float x = curand_uniform(s);
	return x;
}

void DoOneSimulation(SimulationStruct* simulation, int index, char* output, char* fiber1_position) //Wang modified
{
	unsigned long long seed = time(NULL);
	float reflectance[NUM_OF_DETECTOR] = { 0 }; //float reflectance[8] = {0}; //YU-Modified

	MemStruct DeviceMem;
	MemStruct HostMem;
	unsigned int threads_active_total = 1;
	unsigned int i, ii;

	cudaError_t cudastat;

	InitMemStructs(&HostMem, &DeviceMem, simulation, fiber1_position); //Wang modified
	InitDCMem(simulation);

	dim3 dimBlock(NUM_THREADS_PER_BLOCK);	printf("NUM_THREADS_PER_BLOCK\t%d\n", NUM_THREADS_PER_BLOCK);
	dim3 dimGrid(NUM_BLOCKS);				printf("NUM_BLOCKS\t%d\n", NUM_BLOCKS);

	LaunchPhoton_Global << <dimGrid, dimBlock >> >(DeviceMem, seed);
	cudaThreadSynchronize(); //CUDA_SAFE_CALL( cudaThreadSynchronize() ); // Wait for all threads to finish
	cudastat = cudaGetLastError(); // Check if there was an error
	if (cudastat)printf("Error code=%i, %s.\n", cudastat, cudaGetErrorString(cudastat));

	i = 0;

	while (threads_active_total>0)
	{
		i++;
		fiber_initialization(HostMem.f, atof(fiber1_position)); //Wang modified
																//printf("Size of Fibers\t%d\n",sizeof(Fibers));
		cudaMemcpy(DeviceMem.f, HostMem.f, NUM_THREADS * sizeof(Fibers), cudaMemcpyHostToDevice); //malloc sizeof(FIbers) equals to 13*(5*4)

																								  //run the kernel
		seed = time(NULL);
		MCd << <dimGrid, dimBlock >> >(DeviceMem, seed);
		cudaThreadSynchronize(); //CUDA_SAFE_CALL( cudaThreadSynchronize() ); // Wait for all threads to finish
		cudastat = cudaGetLastError(); // Check if there was an error
		if (cudastat)printf("Error code=%i, %s.\n", cudastat, cudaGetErrorString(cudastat));

		// Copy thread_active from device to host, later deleted
		cudaMemcpy(HostMem.thread_active, DeviceMem.thread_active, NUM_THREADS * sizeof(unsigned int), cudaMemcpyDeviceToHost); //CUDA_SAFE_CALL(cudaMemcpy(HostMem.thread_active,DeviceMem.thread_active,NUM_THREADS*sizeof(unsigned int),cudaMemcpyDeviceToHost) );
		threads_active_total = 0;
		for (ii = 0; ii<NUM_THREADS; ii++) threads_active_total += HostMem.thread_active[ii];

		//CUDA_SAFE_CALL(cudaMemcpy(HostMem.num_terminated_photons,DeviceMem.num_terminated_photons,sizeof(unsigned int),cudaMemcpyDeviceToHost) );

		//printf("Run %u, Number of photons terminated %u, Threads active %u\n",i,*HostMem.num_terminated_photons,threads_active_total);

		cudaMemcpy(HostMem.f, DeviceMem.f, NUM_THREADS * sizeof(Fibers), cudaMemcpyDeviceToHost); //CUDA_SAFE_CALL(cudaMemcpy(HostMem.f,DeviceMem.f,NUM_THREADS*sizeof(Fibers),cudaMemcpyDeviceToHost));
		calculate_reflectance(HostMem.f, reflectance);
	}
	//cout << "#" << index << " Simulation done!\n";

	// copy the A_rz and A0_z back to host
	cudaMemcpy(HostMem.A_rz, DeviceMem.A_rz, record_nr * record_nz * sizeof(unsigned long long), cudaMemcpyDeviceToHost);
	cudaMemcpy(HostMem.A0_z, DeviceMem.A0_z, record_nz * sizeof(unsigned long long), cudaMemcpyDeviceToHost);

	output_fiber(simulation, reflectance, output); //Wang modified

	output_A_rz(simulation, HostMem.A_rz, "A_rz.txt"); // output the absorbance
	output_A0_z(simulation, HostMem.A0_z, "A0_z.txt");
	FreeMemStructs(&HostMem, &DeviceMem);
}

void calculate_reflectance(Fibers* f, float *result)
{
	for (int i = 0; i < NUM_THREADS; i++)
	{
		if (NORMAL)
		{
			for (int k = 1; k <= NUM_OF_DETECTOR; k++) {
				result[k - 1] += f[i].data[k];
				/*result[0] += f[i].data[1];
				result[1] += f[i].data[2];
				result[2] += f[i].data[3];
				result[3] += f[i].data[4];*/ //YU-modified
			}
		}
		else
		{   //Wang-modified

			for (int k = 1; k <= NUM_OF_DETECTOR; k++) {
				result[k - 1] += f[i].data[k];
			}

			/*
				// Store the absorbance for grid
				if (p.first_scatter)
				{
					unsigned int index;
					index = min(__float2int_rz(__fdividef(p.z, record_dz)), (int)(record_nz - 1));
					AtomicAddULL(&DeviceMem.A0_z[index], w_temp);
					p.first_scatter = false;
				}
				else
				{
					unsigned int index;
					index = min(__float2int_rz(__fdividef(p.z, record_dz)), (int)(record_nz - 1)) *record_nr + min(__float2int_rz(__fdividef(sqrtf(p.x*p.x + p.y*p.y), record_dr)), (int)record_nr - 1);

					if (index == index_old)
					{
						w += w_temp;
					}
					else// if(w!=0)
					{
						AtomicAddULL(&DeviceMem.A_rz[index_old], w);
						index_old = index;
						w = w_temp;
					}
				}
			*/

		}
	}
}

//Device function to add an unsigned integer to an unsigned long long using CUDA Compute Capability 1.1
__device__ void AtomicAddULL(unsigned long long* address, unsigned int add)
{
	if (atomicAdd((unsigned int*)address, add) + add<add)
		atomicAdd(((unsigned int*)address) + 1, 1u);
}

__global__ void MCd(MemStruct DeviceMem, unsigned long long seed)
{
	//Block index
	int bx = blockIdx.x;

	//Thread index
	int tx = threadIdx.x;

	//First element processed by the block
	int begin = NUM_THREADS_PER_BLOCK * bx;

	float s;	//step length

	unsigned int index, w, index_old;
	index_old = 0;
	w = 0;
	unsigned int w_temp;

	PhotonStruct p = DeviceMem.p[begin + tx];
	Fibers f = DeviceMem.f[begin + tx];

	int new_layer;

	curandState state = DeviceMem.state[begin + tx];
	curand_init(seed, begin + tx, 0, &state);

	//First, make sure the thread (photon) is active
	unsigned int ii = 0;
	if (!DeviceMem.thread_active[begin + tx]) ii = NUMSTEPS_GPU;

	bool k = true;

	for (; ii<NUMSTEPS_GPU; ii++) //this is the main while loop
	{
		if (layers_dc[p.layer].mutr != FLT_MAX)
			s = -__logf(rn_gen(&state))*layers_dc[p.layer].mutr;//sample step length [cm] //HERE AN OPEN_OPEN FUNCTION WOULD BE APPRECIATED
		else
			s = 100.0f;//temporary, say the step in glass is 100 cm.

					   //Check for layer transitions and in case, calculate s
		new_layer = p.layer;
		if (p.z + s*p.dz<layers_dc[p.layer].z_min) {
			new_layer--;
			s = __fdividef(layers_dc[p.layer].z_min - p.z, p.dz);
		} //Check for upwards reflection/transmission & calculate new s
		if (p.z + s*p.dz>layers_dc[p.layer].z_max) {
			new_layer++;
			s = __fdividef(layers_dc[p.layer].z_max - p.z, p.dz);
		} //Check for downward reflection/transmission

		p.x += p.dx*s;
		p.y += p.dy*s;
		p.z += p.dz*s;

		if (p.z>layers_dc[p.layer].z_max)p.z = layers_dc[p.layer].z_max;//needed?
		if (p.z<layers_dc[p.layer].z_min)p.z = layers_dc[p.layer].z_min;//needed?

		if (new_layer != p.layer)
		{
			// set the remaining step length to 0
			s = 0.0f;

			if (Reflect(&p, new_layer, &state) == 0u)//Check for reflection
			{
				if (new_layer == 0)
				{ //Diffuse reflectance
					detect(&p, &f);
					p.weight = 0; // Set the remaining weight to 0, effectively killing the photon
				}
				if (new_layer > *n_layers_dc)
				{	//Transmitted
					p.weight = 0; // Set the remaining weight to 0, effectively killing the photon
				}
			}
		}

		if (s > 0.0f)
		{
			// Drop weight (apparently only when the photon is scattered)
			w_temp = __float2uint_rn(layers_dc[p.layer].mua*layers_dc[p.layer].mutr*__uint2float_rn(p.weight));
			//w_temp = layers_dc[p.layer].mua*layers_dc[p.layer].mutr*p.weight;
			p.weight -= w_temp;
			
			f.absorbed_path[f.absorbed_time][0] = p.x;
			f.absorbed_path[f.absorbed_time][1] = p.y;
			f.absorbed_path[f.absorbed_time][2] = p.z;
			f.absorbed_path[f.absorbed_time][3] = w_temp;
			f.absorbed_time += 1;


			Spin(&p, layers_dc[p.layer].g, &state);
		}

		if (!PhotonSurvive(&p, &state)) //if the photon doesn't survive
		{
			k = false;
			if (atomicAdd(DeviceMem.num_terminated_photons, 1u) < (*num_photons_dc - NUM_THREADS))
			{	// Ok to launch another photon
				LaunchPhoton(&p, &state);//Launch a new photon
			}
			else
			{	// No more photons should be launched. 
				DeviceMem.thread_active[begin + tx] = 0u; // Set thread to inactive
				ii = NUMSTEPS_GPU;				// Exit main loop
			}
		}

	}//end main for loop!
	if (w != 0)
		AtomicAddULL(&DeviceMem.A_rz[index_old], w);

	if (k == true && DeviceMem.thread_active[begin + tx] == 1u)    // photons are not killed after numerous steps
	{
		if (atomicAdd(DeviceMem.num_terminated_photons, 1u) >= (*num_photons_dc - NUM_THREADS))
			DeviceMem.thread_active[begin + tx] = 0u;
	}

	__syncthreads();//necessary?

					//save the state of the MC simulation in global memory before exiting
	DeviceMem.p[begin + tx] = p;	//This one is incoherent!!!
	DeviceMem.f[begin + tx] = f;

}//end MCd

__device__ void LaunchPhoton(PhotonStruct* p, curandState *state)
{
	float rnd_position, rnd_Azimuth, rnd_direction, rnd_rotated;
	float AzimuthAngle;
	float launchPosition;
	float theta_direction;
	float rotated_angle;
	float uxprime, uyprime, uzprime;
	float angle = -ANGLE * PI / 180;


	//rnd_position = rn_gen(state);
	//rnd_Azimuth = rn_gen(state);
	//rnd_direction = rn_gen(state);
	//rnd_rotated = rn_gen(state);
	//AzimuthAngle = 2 * PI * rnd_Azimuth;
	//rotated_angle = 2 * PI * rnd_rotated; //YU-modified // modified to point source

	//float beam_width = 0.0175;  // 175 um, Gaussian beam profile
	//float beam_width = illumination_r;

	// infinite narrow beam for impulse response

	p->x = 0.0;
	p->y = 0.0;
	p->z = 0.0;

	//theta_direction = asin(NAOfSource / n_source)*rnd_direction;
	//p->dz = cos(theta_direction);
	//p->dx = sin(theta_direction) * cos(rotated_angle);
	//p->dy = sin(theta_direction) * sin(rotated_angle);

	//uxprime = cos(angle)*p->dx - sin(angle)*p->dz;
	//uyprime = sin(theta_direction)*sin(rotated_angle);
	//uzprime = sin(angle)*p->dx + cos(angle)*p->dz; // YU-modified

	//p->dx = uxprime, p->dy = uyprime, p->dz = uzprime;
	p->dx = 0.0, p->dy = 0.0, p->dz = 1.0;

	p->layer = 1;
	p->first_scatter = true;

	p->weight = *start_weight_dc; //specular reflection!
}


__global__ void LaunchPhoton_Global(MemStruct DeviceMem, unsigned long long seed)
{
	int bx = blockIdx.x;
	int tx = threadIdx.x;

	//First element processed by the block
	int begin = NUM_THREADS_PER_BLOCK*bx;

	PhotonStruct p;

	curandState state = DeviceMem.state[begin + tx];
	curand_init(seed, 0, 0, &state);

	LaunchPhoton(&p, &state);

	//__syncthreads();//necessary?
	DeviceMem.p[begin + tx] = p;//incoherent!?
}


__device__ void Spin(PhotonStruct* p, float g, curandState *state)
{
	float theta, cost, sint;	// cosine and sine of the 
								// polar deflection angle theta. 
	float cosp, sinp;	// cosine and sine of the 
						// azimuthal angle psi. 
	float temp;
	float tempdir = p->dx;

	//This is more efficient for g!=0 but of course less efficient for g==0
	temp = __fdividef((1.0f - (g)*(g)), (1.0f - (g)+2.0f*(g)*rn_gen(state)));//Should be close close????!!!!!
	cost = __fdividef((1.0f + (g)*(g)-temp*temp), (2.0f*(g)));
	if (g == 0.0f)
		cost = 2.0f*rn_gen(state) - 1.0f;//Should be close close??!!!!!

	sint = sqrtf(1.0f - cost*cost);

	__sincosf(2.0f*PI*rn_gen(state), &sinp, &cosp);// spin psi [0-2*PI)

	temp = sqrtf(1.0f - p->dz*p->dz);

	if (temp == 0.0f) //normal incident.
	{
		p->dx = sint*cosp;
		p->dy = sint*sinp;
		p->dz = copysignf(cost, p->dz*cost);
	}
	else // regular incident.
	{
		p->dx = __fdividef(sint*(p->dx*p->dz*cosp - p->dy*sinp), temp) + p->dx*cost;
		p->dy = __fdividef(sint*(p->dy*p->dz*cosp + tempdir*sinp), temp) + p->dy*cost;
		p->dz = -sint*cosp*temp + p->dz*cost;
	}

	//normalisation seems to be required as we are using floats! Otherwise the small numerical error will accumulate
	temp = rsqrtf(p->dx*p->dx + p->dy*p->dy + p->dz*p->dz);
	p->dx = p->dx*temp;
	p->dy = p->dy*temp;
	p->dz = p->dz*temp;
}// end Spin



__device__ unsigned int Reflect(PhotonStruct* p, int new_layer, curandState *state)
{
	//Calculates whether the photon is reflected (returns 1) or not (returns 0)
	// Reflect() will also update the current photon layer (after transmission) and photon direction (both transmission and reflection)

	float n1 = layers_dc[p->layer].n;
	float n2 = layers_dc[new_layer].n;
	float r;
	float cos_angle_i = fabsf(p->dz);

	if (n1 == n2)//refraction index matching automatic transmission and no direction change
	{
		p->layer = new_layer;
		return 0u;
	}

	if (n1>n2 && n2*n2<n1*n1*(1 - cos_angle_i*cos_angle_i))//total internal reflection, no layer change but z-direction mirroring
	{
		p->dz *= -1.0f;
		return 1u;
	}

	if (cos_angle_i == 1.0f)//normal incident
	{
		r = __fdividef((n1 - n2), (n1 + n2));
		if (rn_gen(state) <= r*r)
		{
			//reflection, no layer change but z-direction mirroring
			p->dz *= -1.0f;
			return 1u;
		}
		else
		{	//transmission, no direction change but layer change
			p->layer = new_layer;
			return 0u;
		}
	}

	//gives almost exactly the same results as the old MCML way of doing the calculation but does it slightly faster
	// save a few multiplications, calculate cos_angle_i^2;
	float e = __fdividef(n1*n1, n2*n2)*(1.0f - cos_angle_i*cos_angle_i); //e is the sin square of the transmission angle
	r = 2 * sqrtf((1.0f - cos_angle_i*cos_angle_i)*(1.0f - e)*e*cos_angle_i*cos_angle_i);//use r as a temporary variable
	e = e + (cos_angle_i*cos_angle_i)*(1.0f - 2.0f*e);//Update the value of e
	r = e*__fdividef((1.0f - e - r), ((1.0f - e + r)*(e + r)));//Calculate r	

	if (rn_gen(state) <= r)
	{
		// Reflection, mirror z-direction!
		p->dz *= -1.0f;
		return 1u;
	}
	else
	{
		// Transmission, update layer and direction
		r = __fdividef(n1, n2);
		e = r*r*(1.0f - cos_angle_i*cos_angle_i); //e is the sin square of the transmission angle
		p->dx *= r;
		p->dy *= r;
		p->dz = copysignf(sqrtf(1 - e), p->dz);
		p->layer = new_layer;
		return 0u;
	}

}

__device__ unsigned int PhotonSurvive(PhotonStruct* p, curandState *state)
{
	//Calculate wether the photon survives (returns 1) or dies (returns 0)

	if (p->weight>WEIGHTI) return 1u; // No roulette needed
	if (p->weight == 0u) return 0u;	// Photon has exited slab, i.e. kill the photon

	if (rn_gen(state) < CHANCE)
	{
		p->weight = __float2uint_rn(__fdividef((float)p->weight, CHANCE));
		//p->weight = __fdividef((float)p->weight,CHANCE);
		return 1u;
	}
	return 0u;
}

__device__ void detect(PhotonStruct* p, Fibers* f)
{
	float angle = ANGLE*PI / 180; //YU-modified
	float critical = asin(f->NA[1] / n_detector); //YU-modified
	float uz_rotated = (p->dx*sin(angle)) + (p->dz*cos(angle)); //YU-modified
	float uz_angle = acos(fabs(uz_rotated)); //YU-modified
	float distance;

	// disable NA consideration //YU-modified
	/*
	distance = sqrt(p->x * p->x + p->y * p->y);

	for (int i = 1; i <= NUM_OF_DETECTOR; i++)
	{
	if ((distance > (f->position[i] - f->radius[i])) && (distance <= (f->position[i] + f->radius[i])))
	{
	//record circular instead of a fiber area
	f->data[i] += p->weight;
	}
	}

	*/

	// NA consideration
	if (uz_angle <= critical)  // successfully detected
	{
		if (NORMAL)
		{
			distance = sqrt(p->x * p->x + p->y * p->y);
			if (distance >= 0.025 && distance <= 0.035)           // SDS = 0.03 cm
				f->data[1] += p->weight / 6.0;
			if (distance >= 0.03 && distance <= 0.05)             // SDS = 0.04 cm
				f->data[2] += p->weight / 16.0;
			if (distance >= 0.05 && distance <= 0.07)             // SDS = 0.06 cm
				f->data[3] += p->weight / 24.0;
			if (distance >= 0.07 && distance <= 0.09)             // SDS = 0.08 cm
				f->data[4] += p->weight / 32.0;
		}
		else
		{
			/*
			for(int i = 1; i <= 3 ; i++)
			{
			if(pow((p->x-f->position[i])*cos(angle),2) + pow(p->y,2) <= f->radius[i]*f->radius[i])
			f->data[i] += p->weight;
			}

			for(int i = 4; i <= 6 ; i++)
			{
			if(pow((p->y-f->position[i]),2) + pow(p->x*cos(angle),2) <= f->radius[i]*f->radius[i])
			f->data[i] += p->weight/2;
			}

			for(int i = 7; i <= 9 ; i++)
			{
			if(pow((p->y-f->position[i]),2) + pow(p->x*cos(angle),2) <= f->radius[i]*f->radius[i])
			f->data[i] += p->weight/2;
			}
			*/
			/*
			distance = sqrt(p->x * p->x + p->y * p->y);
			if (distance >= 0.012 && distance <= 0.032)           // SDS = 0.022 cm
			f->data[1] += p->weight / 8.8;
			if (distance >= 0.035 && distance <= 0.055)             // SDS = 0.045 cm
			f->data[2] += p->weight / 18.0;
			if (distance >= 0.063 && distance <= 0.083)             // SDS = 0.073 cm
			f->data[3] += p->weight / 29.2;

			if (distance >= 0.0115 && distance <= 0.0315)             // SDS = 0.0215 cm
			f->data[4] += p->weight / 8.6;
			if (distance >= 0.031 && distance <= 0.051)             // SDS = 0.041 cm
			f->data[5] += p->weight / 16.4;
			if (distance >= 0.051 && distance <= 0.071)             // SDS = 0.061 cm
			f->data[6] += p->weight / 24.4;
			*/

			distance = sqrt(p->x * p->x + p->y * p->y);
			//all circular
			/*
			for (int i = 1; i <= NUM_OF_DETECTOR; i++)
			{
			if ((distance > (f->position[i] - f->radius[i])) && (distance <= (f->position[i] + f->radius[i])))
			{
			//record circular instead of a fiber area
			f->data[i] += p->weight;
			}
			}
			*/
			//fiber only
			for (int i = 1; i <= NUM_OF_DETECTOR; i++)
			{
				if ((distance >= (f->position[i] - f->radius[i])) && (distance <= (f->position[i] + f->radius[i])))
				{
					float temp;
					temp = (distance*distance + f->position[i] * f->position[i] - f->radius[i] * f->radius[i]) / (2 * distance*f->position[i]);
					// check for rounding error!
					if (temp > 1.0f)
						temp = 1.0f;

					f->data[i] += p->weight  * acos(temp) * RPI;
					f->detected = 1;
				}
			}
		}
	}
	return;
}

int InitDCMem(SimulationStruct* sim)
{
	// Copy num_photons_dc to constant device memory
	cudaMemcpyToSymbol(n_layers_dc, &(sim->n_layers), sizeof(unsigned int));

	// Copy start_weight_dc to constant device memory
	cudaMemcpyToSymbol(start_weight_dc, &(sim->start_weight), sizeof(unsigned int));

	// Copy layer data to constant device memory
	cudaMemcpyToSymbol(layers_dc, sim->layers, (sim->n_layers + 2) * sizeof(LayerStruct));

	// Copy num_photons_dc to constant device memory
	cudaMemcpyToSymbol(num_photons_dc, &(sim->number_of_photons), sizeof(unsigned long long));

	return 0;
}

int InitMemStructs(MemStruct* HostMem, MemStruct* DeviceMem, SimulationStruct* sim, char* fiber1_position) //Wang modified
{
	// Allocate p on the device!!
	cudaMalloc((void**)&DeviceMem->p, NUM_THREADS * sizeof(PhotonStruct));

	// Allocate thread_active on the device and host
	HostMem->thread_active = (unsigned int*)malloc(NUM_THREADS * sizeof(unsigned int));
	if (HostMem->thread_active == NULL) { printf("Error allocating HostMem->thread_active"); exit(1); }
	for (int i = 0; i<NUM_THREADS; i++)HostMem->thread_active[i] = 1u;

	cudaMalloc((void**)&DeviceMem->thread_active, NUM_THREADS * sizeof(unsigned int));
	cudaMemcpy(DeviceMem->thread_active, HostMem->thread_active, NUM_THREADS * sizeof(unsigned int), cudaMemcpyHostToDevice);

	//Allocate num_launched_photons on the device and host
	HostMem->num_terminated_photons = (unsigned long long*) malloc(sizeof(unsigned long long));
	if (HostMem->num_terminated_photons == NULL) { printf("Error allocating HostMem->num_terminated_photons"); exit(1); }
	*HostMem->num_terminated_photons = 0;

	cudaMalloc((void**)&DeviceMem->num_terminated_photons, sizeof(unsigned long long));
	cudaMemcpy(DeviceMem->num_terminated_photons, HostMem->num_terminated_photons, sizeof(unsigned long long), cudaMemcpyHostToDevice);

	//Allocate and initialize fiber f on the device and host
	HostMem->f = (Fibers*)malloc(NUM_THREADS * sizeof(Fibers));
	cudaMalloc((void**)&DeviceMem->f, NUM_THREADS * sizeof(Fibers));
	fiber_initialization(HostMem->f, atof(fiber1_position)); //Wang modified
	cudaMemcpy(DeviceMem->f, HostMem->f, NUM_THREADS * sizeof(Fibers), cudaMemcpyHostToDevice);

	//Allocate states on the device and host
	cudaMalloc((void**)&DeviceMem->state, NUM_THREADS * sizeof(curandState));

	int rz_size;
	rz_size = record_nr*record_nz;

	// Allocate A_rz on host and device
	HostMem->A_rz = (unsigned long long*) malloc(rz_size * sizeof(unsigned long long));
	if (HostMem->A_rz == NULL) { printf("Error allocating HostMem->A_rz"); exit(1); }
	cudaMalloc((void**)&DeviceMem->A_rz, rz_size * sizeof(unsigned long long));
	cudaMemset(DeviceMem->A_rz, 0, rz_size * sizeof(unsigned long long));

	// Allocate A0_z on host and device
	HostMem->A0_z = (unsigned long long*) malloc(record_nz * sizeof(unsigned long long));
	if (HostMem->A0_z == NULL) { printf("Error allocating HostMem->A_rz"); exit(1); }
	cudaMalloc((void**)&DeviceMem->A0_z, record_nz * sizeof(unsigned long long));
	cudaMemset(DeviceMem->A0_z, 0, record_nz * sizeof(unsigned long long));

	return 1;
}

void FreeMemStructs(MemStruct* HostMem, MemStruct* DeviceMem)
{
	free(HostMem->thread_active);
	free(HostMem->num_terminated_photons);
	free(HostMem->f);
	free(HostMem->A_rz);
	free(HostMem->A0_z);

	cudaFree(DeviceMem->thread_active);
	cudaFree(DeviceMem->num_terminated_photons);
	cudaFree(DeviceMem->f);
	cudaFree(DeviceMem->state);
	cudaFree(DeviceMem->A_rz);
	cudaFree(DeviceMem->A0_z);
}