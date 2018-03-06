# CAPE
Cylinder and Plane Extraction from Depth Cameras

## Dependencies

* OpenCV
* Eigen3

## Data

RGBD Sequences are available: [Here](https://drive.google.com/file/d/1v1gNg6JTS3eD_FelXjo9D0yWYaTu7Td6/view?usp=sharing)

Download the .zip file and unzip it into ``./Data``

## Ubuntu Instructions
Tested with Ubuntu 14.04

To compile, inside the directory ``./CAPE`` type:
```
mkdir build
cd build
cmake ../
make
```
To run the executable, type:

```./cape_offline```

## Windows Instructions

Tested configuration: Windows 8.1 with Visual Studio 12

This version includes already a VC11 project.
Just make the necessary changes to link the project with OpenCV and Eigen.
