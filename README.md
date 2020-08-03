# PEN Geant4 software
GEANT4 Simulation to study the PEN optical properties.

It simulates setups that have been mounted at the Max Planck Institute for Physics in Munich in order to study and extract optical paramters of PEN. 

# Compilation
In order to compile the software you need to run within the legend container which contains all the required libraries, e.g Geant4
Then once the software is in your local repository follow these intructions
1) ```venv legeng```
2) create build folder: ```mkdir build``` (if the folder exist you should clean it, ```cd build```, ```rm -rf *```)
3) ```cmake ../```
In case you get errors in this step, you can try: ```cmake -DGeant4_DIR=/opt/geant4/lib64/Geant4-10.3.0 ..```
4) ```make```
This should reate the PEN executable. Then to run in interactive the simulations, just do inside build: ```./PEN```

Select the setup that you want to simulate with:

``` /PEN/det/setDetectorType 1```
Options are:

0 setup with 5 PMTs1 1 PEN sample and no collimator

1 setup for attenuation measurement purposes, consisting in PEN samples with 2 PMTs

2 Same as 1 but with more PEN samples that the user can choose

3 Spectrometer setup. Modified geant4 to use one PMT as spectrometer. It stores the detected photons with its wl information and the traveled distance by the photons before being detected
4 Same as 2 but with collimator

Once you have selected the type of setup that you want to simulate you can chosse the position in which you want to put the collimator together with the trigger setup

```/PEN/det/setCollimatorPositionX 0. mm```

Then you can select the type of material for the target. At the moment it has been implemented: PEN, PVT_structure

```/PEN/det/setTargetMat [material]```

If you select PEN as target material you can also select the light yield with

```/PEN/det/setLY 5000.```

You can also select the roughness of the surfaces with

```/PEN/det/setABS 5.``` 

where the 5. is used as factor to normilize the attenuation curve, which depends of the wl, so you are multiplying by 5. this curve in this case 
Then you can select the number of sampels that you are using in the setup

```/PEN/det/setNTargetSamples 2```

The next step is choose the size fo the target, which depends on the type of detector, in case of setups with 5 PMTs only the thickness has some effect since the the length and width are set by default to 30 mm

```/PEN/det/setTargetThickness 1.7 mm```

In case you are simulating the attenuation setup, you can choose all the dimenssions:

```/PEN/det/setTargetLength 74 mm```

```/PEN/det/setTargetWidth 9 mm```

If you select PEN as material

Then you can select the type of source tou want to simulate with

```/PEN/gun/sourceType 1```

Allowed cases:

  0 - 137Cs source
  
  1 - 207Bi source
  
  2 - 90Sr source
  
  3 - Perpendicular mono-energetic electrons with fixed position - you can change the energy of electrons with 
  
```/PEN/gun/sourceEnergy [number with unit]```

Then you can also select the position of the source with

```/PEN/gun/sourcePositionX 0. mm```

The same applies for z. If you are using the collimator, both position in x,z should be the same, otherwise it won't work properly.

