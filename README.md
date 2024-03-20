# Active Exploration and Reconstruction

This repository is related to an unpublished manuscript.

Active exploration and reconstruction of vascular networks using microrobot swarms demands the cooperation between magnetic actuation platform, imaging devices and microrobot swarms. Here we provide a simple demonstration that can run on a Windows-based environment showing the procedures of decision making, data recording and reconstruction process using recorded top-view and side-view videos. The source code for active exploration and reconstruction are provded for code review as well. More complicated demonstrations demand magnetic actuation of the microrobot swarms and decision-making based on the imaging feedback in real-time. 

## System Requirements
For running the demonstrations:
 - Windows 10
 - MATLAB R2020b
 - LabView 2019 (version 19.0f2 32-bit)
     - Vision acquisition module (version 19.0)
     - Vision development module (version 19.0.0)
     - Mathscript module (version 19.0.0)

For reproducing the experimental results of active exploration using microrobot swarms:
 - Windows 10
 - MATLAB R2020b
 - LabView 2019 (version 19.0f2 32-bit)
     - Vision acquisition module (version 19.0)
     - Vision development module (version 19.0.0)
     - Mathscript module (version 19.0.0)
     - myRIO toolkit (version 7.0)
     - Robotics module (version 19.0.0)
     - Real-time module (version 19.0.0)
 - Hardware platform
     - Dobot CR16 robotic arm
     - Electromagnetic coils with power amplifiers
     - NI myROI microcontroller
     - Top-view and side-view cameras
     - Fe3O4 nanoparticles with SiO2 coating

## Demontration for Active Exploration
Instructions:
1. 

## Reconstruction
Source code for reconstruction process is provided. Detailed descriptions related to the program can be found in our related paper (under review). This program builds the vascular structure using the branch map and reconstruction data recorded from the demonstration in "Demo for Active Exploration\Demo from sample video.vi".
Instructions:
1. 

## Exploration with Data Recording
Source code for real-time mapping, decision making and data recording in active exporation process is provided. The code is in the format of MATLAB function and can be integrated to main program (LabView project in this work) when demanded. Demonstration of this function please refer to "Demo for Active Exploration\Demo from sample video.vi".
