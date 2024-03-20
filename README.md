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

## Sample Code and Demonstrations
### Demontration for Active Exploration
A demonstration showing the image processing, decision making and data recording procedures using the recorded videos are provided in "Demo for Active Exploration\Demo from sample video.vi". 

Instructions:
1. Install the demanded software and modules;
2. Open the demonstration program in "Demo for Active Exploration\Demo from sample video.vi";
3. Run the vi;
4. Click "Start" button to start processing the videos;
5. Results are saved in "Demo for Active Exploration\sample result" in csv format.

### Reconstruction
Source code for reconstruction process is provided in "Reconstruction\Reconstruction.m". Detailed descriptions related to the program can be found in our related paper (under review). This program builds the vascular structure using the branch map and reconstruction data recorded from the demonstration in "Demo for Active Exploration\Demo from sample video.vi".

Instructions:
1. Save the results from demonstration program in "Demo for Active Exploration\sample result" as xlsx format;
2. Open MATLAB and run the program in "Reconstruction\Reconstruction.m";
3. Reconstruction results are shown in figures.

### Exploration with Data Recording
Source code for real-time mapping, decision making and data recording in active exporation process is provided in "Exploration with Data Recording\ActiveExploration_function.m" for code review. The code is in the format of MATLAB function and can be integrated to main program (LabView project in this work) when demanded. Demonstration of this function please refer to "Demo for Active Exploration\Demo from sample video.vi".
