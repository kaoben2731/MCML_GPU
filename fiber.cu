#include "header.h"

void fiber_initialization(Fibers* f, float fiber1_position)
{
	for (int i = 0; i < NUM_THREADS; i++)
	{
		for (int j = 0; j <= NUM_OF_DETECTOR; j++)
			f[i].data[j] = 0;

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


			//f[i].radius[1] = collect_r;					//f[i].radius[1]   = collect_r/2;    //YU-modified         
			//f[i].NA[1]       = NAOfDetector;				
			//f[i].position[1] = 0.0125;                    //first fiber, SDS = 0.03 cm
			//f[i].angle[1]    = ANGLE*PI/180;

			//f[i].radius[2]   = collect_r;          
			//f[i].NA[1]       = NAOfDetector;				
			//f[i].position[2] = 0.0375;                    //second fiber, SDS = 0.04 cm
			//f[i].angle[2]    = ANGLE*PI/180;		
			//
			//f[i].radius[3]   = collect_r;             
			//f[i].NA[3]       = NAOfDetector;				
			//f[i].position[3] = 0.0625;                    //third fiber, SDS = 0.06 cm
			//f[i].angle[3]    = ANGLE*PI/180;

			//f[i].radius[4]   = collect_r;             
			//f[i].NA[4]       = NAOfDetector;				
			//f[i].position[4] = 0.0875;                    //fourth fiber, SDS = 0.08 cm
			//f[i].angle[4]    = ANGLE*PI/180;

			//f[i].radius[5] = collect_r;
			//f[i].NA[5] = NAOfDetector;
			//f[i].position[5] = 0.15;                    //fourth fiber, SDS = 0.08 cm
			//f[i].angle[5] = ANGLE*PI / 180;

			//f[i].radius[6] = collect_r;
			//f[i].NA[6] = NAOfDetector;
			//f[i].position[6] = 0.20;                    //fourth fiber, SDS = 0.08 cm
			//f[i].angle[6] = ANGLE*PI / 180;
		}
		else
		{
			/*
			for (int k = 1; k <= NUM_OF_DETECTOR; k++) {
			f[i].radius[k] = collect_r;
			f[i].NA[k] = NAOfDetector;
			f[i].position[k] = 0.74 + 0.02*(k-1);
			f[i].angle[k] = ANGLE*PI / 180;
			}
			*/
			/*
			f[i].radius[1] = collect_r;      //fiber 1
			f[i].NA[1] = NAOfDetector;
			f[i].position[1] = 0.022;
			f[i].angle[1] = ANGLE*PI / 180;

			f[i].radius[2] = collect_r;      //fuber 2
			f[i].NA[2] = NAOfDetector;
			f[i].position[2] = 0.045;
			f[i].angle[2] = ANGLE*PI / 180;

			f[i].radius[3] = collect_r;      //fiber 3
			f[i].NA[3] = NAOfDetector;
			f[i].position[3] = 0.073;
			f[i].angle[3] = ANGLE*PI / 180;
			*/
			//TiO2 test
			// muscle fiber
			/*
			f[i].radius[1] = collect_r;      //fiber 1
			f[i].NA[1] = NAOfDetector;
			f[i].position[1] = 0.203;
			f[i].angle[1] = ANGLE*PI / 180;

			f[i].radius[2] = collect_r;      //fuber 2
			f[i].NA[2] = NAOfDetector;
			f[i].position[2] = 0.260;
			f[i].angle[2] = ANGLE*PI / 180;

			f[i].radius[3] = collect_r;      //fiber 3
			f[i].NA[3] = NAOfDetector;
			f[i].position[3] = 0.302;
			f[i].angle[3] = ANGLE*PI / 180;

			f[i].radius[4] = collect_r;      //fiber 4
			f[i].NA[4] = NAOfDetector;
			f[i].position[4] = 0.363;
			f[i].angle[4] = ANGLE*PI / 180;
			*/
			/*
			for (int k = 1; k <= NUM_OF_DETECTOR; k++) {
			f[i].radius[k] = collect_r;
			f[i].NA[k] = NAOfDetector;
			f[i].position[k] = fiber1_position + 0.2*(k - 1);
			f[i].angle[k] = ANGLE*PI / 180;
			}
			*/
			//for skin - collagen
			/*
			for (int k = 1; k <= NUM_OF_DETECTOR; k++) {
			f[i].radius[k] = collect_r;
			f[i].NA[k] = NAOfDetector;
			f[i].position[k] = fiber1_position + 0.1*(k - 1);
			f[i].angle[k] = ANGLE*PI / 180;
			}
			*/
			/*
			//for skin - melanin
			f[i].radius[1] = collect_r;      //fiber 1
			f[i].NA[1] = NAOfDetector;
			f[i].position[1] = 0.05;
			f[i].angle[1] = ANGLE*PI / 180;

			f[i].radius[2] = collect_r;      //fiber 2
			f[i].NA[2] = NAOfDetector;
			f[i].position[2] = 0.10;
			f[i].angle[2] = ANGLE*PI / 180;

			f[i].radius[3] = collect_r;      //fiber 3
			f[i].NA[3] = NAOfDetector;
			f[i].position[3] = 0.15;
			f[i].angle[3] = ANGLE*PI / 180;

			for (int k = 4; k <= NUM_OF_DETECTOR; k++) {
			f[i].radius[k] = collect_r;
			f[i].NA[k] = NAOfDetector;
			f[i].position[k] = fiber1_position + 0.1*(k - 4);
			f[i].angle[k] = ANGLE*PI / 180;
			}
			*/
			//IJV holder
			/*
			for (int k = 1; k <= NUM_OF_DETECTOR; k++) {
			f[i].radius[k] = collect_r;
			f[i].NA[k] = NAOfDetector;
			f[i].position[k] = fiber1_position + 0.2*(k - 1);
			f[i].angle[k] = ANGLE*PI / 180;
			}
			*/

			/*
			for (int k = 1; k <= NUM_OF_DETECTOR; k++) {
			f[i].radius[k] = collect_r;
			f[i].NA[k] = NAOfDetector;
			f[i].position[k] = fiber1_position + 0.05*(k - 1);
			f[i].angle[k] = ANGLE*PI / 180;
			}
			*/

			// for Apacer skin exp
			/*
			f[i].radius[1] = collect_r;      //fiber 1
			f[i].NA[1] = NAOfDetector;
			f[i].position[1] = 0.0845;
			f[i].angle[1] = ANGLE*PI / 180;

			f[i].radius[2] = collect_r;      //fuber 2
			f[i].NA[2] = NAOfDetector;
			f[i].position[2] = 0.1455;
			f[i].angle[2] = ANGLE*PI / 180;

			f[i].radius[3] = collect_r;      //fiber 3
			f[i].NA[3] = NAOfDetector;
			f[i].position[3] = 0.1875;
			f[i].angle[3] = ANGLE*PI / 180;

			f[i].radius[4] = collect_r;      //fiber 4
			f[i].NA[4] = NAOfDetector;
			f[i].position[4] = 0.2445;
			f[i].angle[4] = ANGLE*PI / 180;

			f[i].radius[5] = collect_r;      //fiber 5
			f[i].NA[5] = NAOfDetector;
			f[i].position[5] = 0.3845;
			f[i].angle[5] = ANGLE*PI / 180;

			f[i].radius[6] = collect_r;      //fiber 6
			f[i].NA[6] = NAOfDetector;
			f[i].position[6] = 0.4455;
			f[i].angle[6] = ANGLE*PI / 180;

			f[i].radius[7] = collect_r;      //fiber 7
			f[i].NA[7] = NAOfDetector;
			f[i].position[7] = 0.4875;
			f[i].angle[7] = ANGLE*PI / 180;

			f[i].radius[8] = collect_r;      //fiber 8
			f[i].NA[8] = NAOfDetector;
			f[i].position[8] = 0.5445;
			f[i].angle[8] = ANGLE*PI / 180;

			f[i].radius[9] = collect_r;      //fiber 9
			f[i].NA[9] = NAOfDetector;
			f[i].position[9] = 0.6845;
			f[i].angle[9] = ANGLE*PI / 180;

			f[i].radius[10] = collect_r;      //fiber 10
			f[i].NA[10] = NAOfDetector;
			f[i].position[10] = 0.7455;
			f[i].angle[10] = ANGLE*PI / 180;

			f[i].radius[11] = collect_r;      //fiber 11
			f[i].NA[11] = NAOfDetector;
			f[i].position[11] = 0.7875;
			f[i].angle[11] = ANGLE*PI / 180;

			f[i].radius[12] = collect_r;      //fiber 12
			f[i].NA[12] = NAOfDetector;
			f[i].position[12] = 0.8445;
			f[i].angle[12] = ANGLE*PI / 180;
			*/
			// for New muscle probe
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

		}

	}
}