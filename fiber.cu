#include "header.h"

void fiber_initialization(Fibers* f, int num_of_threads)
{
	for (int i = 0; i < num_of_threads; i++)
	{
		for (int j = 0; j < SDS_detected_temp_size; j++) {
			f[i].data[j] = 0;
			f[i].detected_SDS_number[j] = 0;
		}
		f[i].detected_photon_counter = 0;
	}
}

void fiber_initialization_replay(Fibers_Replay* f_r, SimulationStruct* sim, int num_of_threads)
{
	for (int i = 0; i < num_of_threads; i++)
	{
		f_r[i].have_detected = false;
		f_r[i].data = 0;
		f_r[i].scatter_event = 0;
		f_r[i].detected_SDS_number = 0;
		for (int l = 0; l < sim->num_layers; l++)
		{
			f_r[i].layer_pathlength[l] = 0;
		}
	}
}