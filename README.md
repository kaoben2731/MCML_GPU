# MCML_GPU
the MCML with cuda acceleration which modified by Benjamin Kao and have many purpose


---

## Contents

* [Prepare](#prepare)
* [Execute the Program](#exec)
* [Input Files](#input_files)
    * [sim_set.json](#simset_json)
    * [input.txt](#input_txt)
* [Output Files](#output_files)
* [output.txt](#output_txt)
* [Summary.json](#summary_json)
* [pathlength_SDS_#.txt or pathlength_SDS_#.bin](#PLSDS_txt)
* [A_rz_SDS_#.txt](#arz_txt)
* [A0_z_SDS_#.txt](#a0z_txt)
* [average_pathlength.txt](#avg_PL_txt)
* [Update History](#update_history)
* [Reference](#reference)

---

<h2 id="prepare">Prepare</h2>

1. Install Nvidia Cuda driver on the computer.
2. type `make` to make the program

---

<h2 id="exec">Execute the Program</h2>

Use this command to run this program:    
`./MCML_GPU sim_set.json input.txt output.txt <option(s)>`
* [`sim_set.json`](#simset_json)
A .json file for setting parameters.
* [`input.txt`](#input_txt)
The input file to set optical parameters for each layer.
* [`output.txt`](#output_txt)
The output filename for reflectance.
* options
    * `-h`
Print the helping information.
    * `-R`
Replay the detected photons after first simulation, to get the pathlength in each layer or absorbance matrix.
	* `-A`
Output the absorbance array for each detector in [A_rz_SDS_#.txt](#arz_txt) and [A0_z_SDS_#.txt](#a0z_txt).
    * `-P`
Output the pathlength for each photon in [pathlength_SDS_#.txt](#PLSDS_txt), otherwise output the calculated average pathlength in [average_pathlength.txt](#avg_PL_txt).
    * `-AP`
Calaulate and output the average pathlength.
    * `-B`
Output the pathlength file in binary format for faster speed.
    * `-G #`
Select which GPU to use.
    * `-LS`
List the GPUs in the computer.


---

<h2 id="input_files">Input Files</h2>

<h3 id="simset_json">sim_set.json</h3>

Use a .json file to set the parameter of simulations
* number_simulation:
How many simulations to run, usually set to 1.
* number_photons:
How many photons in one simulaiton
* number_layers:
How many layer the tissue is
* detector_reflectance:
The reflectance rate of the detector.  If no reflectance, set to 0.
* upper_n:
The refractive index for upper layer outside the tissue
* lower_n:
The refractive index for lower layer outside the tissue
* source_probe_oblique:
The source probe is oblique or not.  If this is 1, then the simulation will use fiber detection mode; otherwise, use ring detection mode.
* detector_probe_oblique:
The detector probes are oblique or not
* probes:
The parameters for the fibers
* source:
Parameters for source fiber
    * NA:
Numerical aperture of source fiber
    * radius:
The radius of source fiber, in cm
    * angle:
The angle of source fiber, in degree, toward +x direction
* num_SDS:
Number of detector fibers
* detectors:
Parameters for detector fibers, shold be an array of the same size as "num_SDS"
    * pos:
The position (distance of the center of source fiber to the center of this detector fiber), in cm
    * NA:
Numerical aperture of this detector fiber
    * radius:
The radius of this detector fiber, in cm
    * angle:
The angle of this detector fiber, in degree, toward -x direction
		
* Example:
```
{
  "number_simulation": 1,
  "number_photons": 100000000,
  "number_layers": 5,
  "detector_reflectance": 0.0,
  "upper_n": 1.457,
  "lower_n": 1.457,
  "source_probe_oblique": 1,
  "detector_probe_oblique": 1,
  "probes": {
    "source": {
      "NA": 0.37,
      "radius": 0.075,
      "angle": 0
    },
    "num_SDS": 3,
    "detectors": [
      {
        "pos": 0.8,
        "NA": 0.13,
        "radius": 0.001,
        "angle": 0
      },
      {
        "pos": 1.5,
        "NA": 0.12,
        "radius": 0.02,
        "angle": 0
      },
      {
        "pos": 3.0,
        "NA": 0.12,
        "radius": 0.02,
        "angle": 0
      }
    ]
  }
}
```

---

<h3 id="input_txt">input.txt</h3>

Set the optical parameters for each layer.

* Rules:
    * Each line for one simulation.
    * Arrange in "height(cm) mu_a(1/cm) mu_s(1/cm) n g" for one layer.
    * The last layer should have no height.
    * 5 parameters for 1 layer, so n layer should be 5n-1 columns.
    * The number of layers should be the same as "number_layers" in `sim_set.json`
* Example:
`0.27	0.3201	35.2861	1.4	0	0.65	0.2954	30.5833	1.4	0	0.2747	130.6813	1.4	0`
This is a 3 layer tissue, parameters for each layer are setting as below:

| layer    | height   | mu_a     | mu_s     | n        | g        |
| -------- | -------- | -------- | -------- | -------- | -------- |
| 1        | 0.27     | 0.3201   | 35.2861  | 1.4      | 0        |
| 2        | 0.65     | 0.2954   | 30.5833  | 1.4      | 0        |
| 3        | N/A      | 0.2747   | 130.6813 | 1.4      | 0        |


---

<h2 id="output_files">Output Files</h2>

<h3 id="output_txt">output.txt</h3>

The reflectance collected by each SDS.  One column is one SDS.

---

<h3 id="summary_json">Summary.json</h3>

Store the summary of this simulaiton.
* num_photon:
Number of total simulated photon.
* sim_time:
The time cost to do the (first) simulate, in secs
* replay_time:
The time cost to do the second simulate and output the pathlength array, in secs
* sim_speed(photons/s):
The average speed of simulation.
* sim_GPU:
The name of GPU.
* each_photon_weight:
The initial weight of each photon at launch, when calculate WMC, the **[photon detected weight](#PLSDS_txt)** should be divided by this number.
* number_layers:
The number of tissue layers.
* detect_mode:
Use ring detector to collect photons, or using little fiber detector.
* num_SDS:
The number of derector fibers.
* SDS_detected_number:
The number of photon detected by each detector.  The number of elements is the same as **num_SDS** in [sim_set.json](#simset_json).

---

<h3 id="PLSDS_txt">pathlength_SDS_#.txt or pathlength_SDS_#.bin</h3>

The pathlength information for each photon collected by SDS #, can be used to perform White MC.
* Rules:
    * Each row is one photon.
    * The first column is the photon's weight when it been detected.
    * The second column to the second-last columns are the photon's pathlength (cm) in each tissue layer.
    * The last column is how many times did the photon scattered.
* Example:
```
7155930	0.164281	2.99355	4.65341	570
6295.91	0.0202396	0.570237	0.312309	3036
827985	0.108664	2.41008	4.92397	1120
325157	0.225788	1.57377	4.2767	873
17125.5	0.0720376	1.84378	2.97216	2504
27530.4	0.106237	1.10873	1.53136	2460
```
Those are 6 detected photon with properties below:

| weight   | PL in L1  | PL in L2 | PL in L3 | scatter  |
| -------- | --------  | -------- | -------- | -------- |
| 7155930  | 0.164281  | 2.99355  | 4.65341  | 570      |
| 6295.91  | 0.0202396 | 0.570237 | 0.312309 | 3036     |
| 827985   | 0.108664  | 2.41008  | 4.92397  | 1120     |
| 325157   | 0.225788  | 1.57377  | 4.2767   | 873      |
| 17125.5  | 0.0720376 | 1.84378  | 2.97216  | 2504     |
| 27530.4  | 0.106237  | 1.10873  | 1.53136  | 2460     |

---

<h3 id="arz_txt">A_rz_SDS_#.txt</h3>

* The absorbance array (except the first scatter event) for SDS #.
* In the format of a 2-D array, with rows are z-directional and columns are r-directional.

---

<h3 id="a0z_txt">A0_z_SDS_#.txt</h3>

* The absorbance array the first scatter event along z-axis for SDS #.
* In the format of a 1-D array, with rows are z-directional.

---

<h3 id="avg_PL_txt">average_pathlength.txt</h3>

* The average pathlength (cm) for each detector in each layer.
* Start from the PL in 1st layer for 1st detector, in 2nd layer for 1st detector......
* If there are 4 detectors and 5 layers of tissue, than there will be 20 columns.

---

<h2 id="update_history">Update History</h2>

---

### K1.01

* updata: 2019/10/12
* Add the `-A` option to control output absorbance array.
* Change the original `-A` (average pathlength) option to `-AP`.

### K2.01

* updata: 2019/12/03
* Add source_probe_oblique, detector_probe_oblique options in the setting file.
* Add the oblique mode in simulation, which will use fiber detection mode.
* Read and output GPU information.
* Set the number of blocks and threads according to the GPU information.

### K2.02
* updata: 2020/02/22
* Add binary output mode (`-B`) for absorbance array (`-A`).

### K2.03
* updata: 2020/02/23
* Implement source NA and radius setting.

### K3.01
* updata: 2020/02/23
* Can choose which GPU to run the simulation (`-G #`)
* Can list the GPUs in the computer (`-LS`)

### K3.02
* updata: 2020/06/13
* Little update in decide the scatter angle of the photon to improve speed.

---

<h2 id="reference">Reference</h2>

### The original CUDAMCML

[code link](http://www.atomic.physics.lu.se/biophotonics/research/monte-carlo-simulations/gpu-monte-carlo/)

Erik Alerstam, Tomas Svensson, and Stefan Andersson-Engels. "Parallel
Computing with Graphics Processing Units for High-Speed Monte Carlo
Simulation of Photon Migration." Journal of Biomedical Optics 13(6), 060504 (2008)


---

<h2 id="todo">TODO</h2>

---

* Check Why I **can't** get the same reflectance while the musp is the same!
* Setting of detector n different from outer n
* Change the hight setting of buttom layer to true Inf
* NA and radius setting for oblique probe
* Improve the oblique probe performance on linux
* set n_dr, n_dz in setup file
* add multiple GPU
* try to improve performance by increasing the num_step
* add a preview function to draw the settings.