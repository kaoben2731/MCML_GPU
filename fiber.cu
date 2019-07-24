#include "header.h"

void fiber_initialization(Fibers* f, SimulationStruct* sim)
{
	for (int i = 0; i < NUM_THREADS; i++)
	{
		for (int j = 0; j < SDS_detected_temp_size; j++) {
			f[i].data[j] = 0;
			f[i].detected_SDS_number[j] = 0;
		}

		f[i].detected_photon_counter = 0;
		
		f[i].radius = new float[sim->num_detector + 1];
		f[i].NA = new float[sim->num_detector + 1];
		f[i].angle = new float[sim->num_detector + 1];
		f[i].position = new float[sim->num_detector + 1];

		for (int d = 0; d <= sim->num_detector; d++) {
			f[i].radius[d] = sim->detInfo[d].raduis;
			f[i].NA[d] = sim->detInfo[d].NA;
			f[i].angle[d] = ANGLE*PI / 180;
			if (d == 0) {
				f[i].position[d] = 0.0;
			}
			else {
				f[i].position[d] = sim->detInfo[d].position;
			}
		}

		// for skin
		//f[i].radius[1] = collect_r;      //fiber 1   // turn 200um into 10um
		//f[i].NA[1] = NAOfDetector;
		//f[i].position[1] = 0.0215;
		//f[i].angle[1] = ANGLE*PI / 180;

		//f[i].radius[2] = collect_r;      //fuber 2
		//f[i].NA[2] = NAOfDetector;
		//f[i].position[2] = 0.022;
		//f[i].angle[2] = ANGLE*PI / 180;

		//f[i].radius[3] = collect_r;      //fiber 3
		//f[i].NA[3] = NAOfDetector;
		//f[i].position[3] = 0.041;
		//f[i].angle[3] = ANGLE*PI / 180;

		//f[i].radius[4] = collect_r;      //fiber 4
		//f[i].NA[4] = NAOfDetector;
		//f[i].position[4] = 0.045;
		//f[i].angle[4] = ANGLE*PI / 180;

		//f[i].radius[5] = collect_r;      //fiber 5
		//f[i].NA[5] = NAOfDetector;
		//f[i].position[5] = 0.061;
		//f[i].angle[5] = ANGLE*PI / 180;

		//f[i].radius[6] = collect_r;      //fiber 6
		//f[i].NA[6] = NAOfDetector;
		//f[i].position[6] = 0.073;
		//f[i].angle[6] = ANGLE*PI / 180;

		// for shorter muscle
		//f[i].radius[1] = collect_r;      //fiber 1   // turn 200um into 10um
		//f[i].NA[1] = NAOfDetector;
		//f[i].position[1] = 0.07;
		//f[i].angle[1] = ANGLE*PI / 180;

		//f[i].radius[2] = collect_r;      //fuber 2
		//f[i].NA[2] = NAOfDetector;
		//f[i].position[2] = 0.1;
		//f[i].angle[2] = ANGLE*PI / 180;

		//f[i].radius[3] = collect_r;      //fiber 3
		//f[i].NA[3] = NAOfDetector;
		//f[i].position[3] = 0.12;
		//f[i].angle[3] = ANGLE*PI / 180;

		//f[i].radius[4] = collect_r;      //fiber 4
		//f[i].NA[4] = NAOfDetector;
		//f[i].position[4] = 0.16;
		//f[i].angle[4] = ANGLE*PI / 180;

		//f[i].radius[5] = collect_r;      //fiber 5
		//f[i].NA[5] = NAOfDetector;
		//f[i].position[5] = 0.2;
		//f[i].angle[5] = ANGLE*PI / 180;

		//f[i].radius[6] = collect_r;      //fiber 6
		//f[i].NA[6] = NAOfDetector;
		//f[i].position[6] = 0.25;
		//f[i].angle[6] = ANGLE*PI / 180;

		/*
		// for muscle
		f[i].radius[1] = collect_r;      //fiber 1   // turn 200um into 10um
		f[i].NA[1] = NAOfDetector;
		f[i].position[1] = 0.314;
		f[i].angle[1] = ANGLE*PI / 180;

		f[i].radius[2] = collect_r;      //fuber 2
		f[i].NA[2] = NAOfDetector;
		f[i].position[2] = 0.401;
		f[i].angle[2] = ANGLE*PI / 180;

		f[i].radius[3] = collect_r;      //fiber 3
		f[i].NA[3] = NAOfDetector;
		f[i].position[3] = 0.444;
		f[i].angle[3] = ANGLE*PI / 180;

		f[i].radius[4] = collect_r;      //fiber 4
		f[i].NA[4] = NAOfDetector;
		f[i].position[4] = 0.535;
		f[i].angle[4] = ANGLE*PI / 180;

		f[i].radius[5] = collect_r;      //fiber 5
		f[i].NA[5] = NAOfDetector;
		f[i].position[5] = 0.626;
		f[i].angle[5] = ANGLE*PI / 180;
		*/

		/*
		// for longer muscle
		f[i].radius[1] = collect_r;      //fiber 1   // turn 200um into 10um
		f[i].NA[1] = NAOfDetector;
		f[i].position[1] = 0.75;
		f[i].angle[1] = ANGLE*PI / 180;

		f[i].radius[2] = collect_r;      //fuber 2
		f[i].NA[2] = NAOfDetector;
		f[i].position[2] = 0.85;
		f[i].angle[2] = ANGLE*PI / 180;

		f[i].radius[3] = collect_r;      //fiber 3
		f[i].NA[3] = NAOfDetector;
		f[i].position[3] = 1;
		f[i].angle[3] = ANGLE*PI / 180;

		f[i].radius[4] = collect_r;      //fiber 4
		f[i].NA[4] = NAOfDetector;
		f[i].position[4] = 1.2;
		f[i].angle[4] = ANGLE*PI / 180;

		f[i].radius[5] = collect_r;      //fiber 5
		f[i].NA[5] = NAOfDetector;
		f[i].position[5] = 1.35;
		f[i].angle[5] = ANGLE*PI / 180;

		f[i].radius[6] = collect_r;      //fiber 6
		f[i].NA[6] = NAOfDetector;
		f[i].position[6] = 1.5;
		f[i].angle[6] = ANGLE*PI / 180;
		*/

		// for UT-Austin's system also our fNIRS system
		//f[i].radius[1] = collect_r/20;      //fiber 1   // turn 200um into 10um
		//f[i].NA[1] = NAOfDetector;
		//f[i].position[1] = 0.8;
		//f[i].angle[1] = ANGLE*PI / 180;

		//f[i].radius[2] = collect_r;      //fuber 2
		//f[i].NA[2] = NAOfDetector;
		//f[i].position[2] = 1.5;
		//f[i].angle[2] = ANGLE*PI / 180;

		//f[i].radius[3] = collect_r;      //fiber 3
		//f[i].NA[3] = NAOfDetector;
		//f[i].position[3] = 3.0;
		//f[i].angle[3] = ANGLE*PI / 180;

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

void fiber_initialization_replay(Fibers_Replay* f_r, SimulationStruct* sim)
{
	for (int i = 0; i < NUM_THREADS; i++)
	{
		f_r[i].have_detected = false;
		f_r[i].data = 0;
		f_r[i].scatter_event = 0;
		f_r[i].detected_SDS_number = 0;
		f_r[i].layer_pathlength = new float[sim->num_layers];
		for (int l = 0; l < sim->num_layers; l++)
		{
			f_r[i].layer_pathlength[l] = 0;
		}
	}
}