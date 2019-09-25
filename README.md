# MCML_GPU
the MCML with cuda acceleration which modified by Benjamin and have many purpose


---

## execute the program
./MCML_GPU sim_set.json input.txt output.txt <option(s)>
### sim_set.json
a .json file for setting parameters
### input.txt
the input file to set optical parameters for each layer
### output.txt
the output filename for reflectance
### options
#### -h
Print the helping information
#### -R
Replay the detected photons after first simulation, to get the pathlength in each layer or absorbance matrix
#### -P
Output the pathlength for each photon, otherwise output the calculated average pathlength
#### -A
Calaulate and output the average pathlength
#### -B
Output the pathlength file in binary format


---

## sim_set.json
use a .json file to set the parameter of simulations
### number_simulation
how many simulations to run, usually set to 1.
### number_photons
how many photons in one simulaiton
### number_layers
how many layer the tissue is
### detector_reflectance
the reflectance rate of the detector.  If no reflectance, set to 0.
### upper_n
the refractive index for upper layer outside the tissue
### lower_n
the refractive index for lower layer outside the tissue
### probes
the parameters for the fibers
#### source
parameters for source fiber
##### NA
numerical aperture of source fiber
##### radius
the radius of source fiber, in cm
#### num_SDS
number of detector fibers
#### detectors
parameters for detector fibers, shold be an array of the same size as "num_SDS"
##### pos
the position (distance of the center of source fiber to the center of this detector fiber), in cm
##### NA
numerical aperture of this detector fiber
##### radius
the radius of this detector fiber, in cm

---

## input.txt
Set the optical parameters for each layer.

### Rules:
* Each line for one simulation.
* Arrange in "height(cm) mu_a(1/cm) mu_s(1/cm) n g" for one layer.
* The last layer should have no height.
* 5 parameters for 1 layer, so n layer should be 5n-1 columns.
* The number of layers should be the same as "number_layers" in 
### Example
`0.27	0.3201	35.2861	1.4	0	0.65	0.2954	30.5833	1.4	0	0.2747	130.6813	1.4	0`
This is a 3 layer tissue, parameters for each layer are setting as below:

| layer    | height   | mu_a     | mu_s     | n        | g        |
| -------- | -------- | -------- | -------- | -------- | -------- |
| 1        | 0.27     | 0.3201   | 35.2861  | 1.4      | 0        |
| 2        | 0.65     | 0.2954   | 30.5833  | 1.4      | 0        |
| 3        | N/A      | 0.2747   | 130.6813 | 1.4      | 0        |


---

## output files
