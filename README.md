# CAPE
Cylinder and Plane Extraction from Depth Cameras

Implementation of the method proposed in:  
Fast Cylinder and Plane Extraction from Depth Cameras for Visual Odometry,   
https://arxiv.org/abs/1803.02380

Note: The parameters are fine-tuned for detecting large surfaces with Kinect 1 and Structure sensor for VO. For other applications, these may need to be modified.

## Dependencies

* OpenCV
* Eigen3

## Data

RGBD Sequences w/ cylinders are available for testing: [Here](https://drive.google.com/drive/folders/1CaVVLF7AQUlsOwFWrx-Fm7zB6wueQBE3?usp=sharing)

Download a .zip file and unzip it into ``./Data``

## Ubuntu Instructions
Tested with Ubuntu 14.04 & 16.04

To compile, inside the directory ``./CAPE`` which contains the CmakeLists.txt, type:
```
mkdir build
cd build
cmake ..
make
```
To run the executable, type:

```./cape_offline 20 <sequence_name>```

where the first argument is the cell size in pixels (recommended 20)
and <sequence_name> is a folder:
``./Data/<sequence_name>``
that contains the image and calibration files

## Windows Instructions

Tested configuration: Windows 8.1 with Visual Studio 12

This version includes already a VC11 project.
Just make the necessary changes to link the project with OpenCV and Eigen.
