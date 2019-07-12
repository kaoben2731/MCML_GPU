#include "header.h"

void fiber_initialization(Fibers* f, float fiber1_position)
{
	for (int i = 0; i < NUM_THREADS; i++)
	{
		for (int j = 0; j < SDS_detected_temp_size; j++) {
			f[i].data[j] = 0;
			f[i].detected_SDS_number[j] = 0;
		}

		f[i].detected_photon_counter = 0;
		
		f[i].radius[0] = illumination_r;          // source fiber			
		f[i].NA[0] = NAOfSource;
		f[i].angle[0] = ANGLE*PI / 180;
		f[i].position[0] = 0.0;

		if (NORMAL)
		{
			for (int k = 1; k <= NUM_OF_DETECTOR; k++) {
				f[i].radius[k] = collect_r;
				f[i].NA[k] = NAOfDetector;
				f[i].position[k] = collect_r + 2.0 * (collect_r) * (float)(k - 1);
				f[i].angle[k] = ANGLE*PI / 180;
				//if(i==0 && k<= NUM_OF_DETECTOR)printf("p[%d][%d]\t%f\n",i,k,f[i].position[k]);
			}
		}
		else
		{
			// for UT-Austin's system also our fNIRS system
			f[i].radius[1] = collect_r/20;      //fiber 1   // turn 200um into 10um
			f[i].NA[1] = NAOfDetector;
			f[i].position[1] = 0.8;
			f[i].angle[1] = ANGLE*PI / 180;

			f[i].radius[2] = collect_r;      //fuber 2
			f[i].NA[2] = NAOfDetector;
			f[i].position[2] = 1.5;
			f[i].angle[2] = ANGLE*PI / 180;

			f[i].radius[3] = collect_r;      //fiber 3
			f[i].NA[3] = NAOfDetector;
			f[i].position[3] = 3.0;
			f[i].angle[3] = ANGLE*PI / 180;

			// for our old NIRS system
			/*
			f[i].radius[1] = collect_r;      //fiber 1
			f[i].NA[1] = NAOfDetector;
			f[i].position[1] = 1.45;
			f[i].angle[1] = ANGLE*PI / 180;

			f[i].radius[2] = collect_r;      //fuber 2
			f[i].NA[2] = NAOfDetector;
			f[i].position[2] = 2.05;
			f[i].angle[2] = ANGLE*PI / 180;

			f[i].radius[3] = collect_r;      //fiber 3
			f[i].NA[3] = NAOfDetector;
			f[i].position[3] = 2.9;
			f[i].angle[3] = ANGLE*PI / 180;

			f[i].radius[4] = collect_r;      //fiber 4
			f[i].NA[4] = NAOfDetector;
			f[i].position[4] = 3.24;
			f[i].angle[4] = ANGLE*PI / 180;

			f[i].radius[5] = collect_r;      //fiber 5
			f[i].NA[5] = NAOfDetector;
			f[i].position[5] = 4.1;
			f[i].angle[5] = ANGLE*PI / 180;
			*/

		}

	}
}

void fiber_initialization_replay(Fibers_Replay* f_r)
{
	for (int i = 0; i < NUM_THREADS; i++)
	{
		f_r[i].have_detected = false;
		f_r[i].data = 0;
		f_r[i].scatter_event = 0;
		f_r[i].detected_SDS_number = 0;
		for (int l = 0; l < NUM_LAYER; l++)
		{
			f_r[i].layer_pathlength[l] = 0;
		}
	}
}