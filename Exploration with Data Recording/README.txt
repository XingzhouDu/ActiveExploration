Copyright © 2023 DU Xingzhou.
Licence: Apache license 2.0.

Source code for real-time mapping, decision making and data recording in active exporation process. The code is in the format of MATLAB function and can be integrated to main program of magnetic actuation system (LabView project in this work) when demanded. The input and output variables are described as follows.

Input:
  imgroi2_top - binarized ROI image on top view with a resolution of 256*256, in which the value of each pixel is 0 or 255;
  img_grey_side - greyscale ROI image on side view with a resolution of 128*400, in which the value of each pixel is 0 to 255;
  img_grey_side_ini - greyscale background side-view image with a resolution of 128*400;
  len_x, len_y - consts showing the dimension of top view (in millimeter) in x and y directions, respectively;
  ROI_len - length of the dynamic ROI on top view in pixels;
  ROI_hei - height of the dynamic ROI on side view in pixels;
  xyROI_p - a 1*4 vector recording the dynamic ROI on top view in last cycle, in the format of [x_upperleft, y_upperleft, x_lowerright, y_lowerright], initialized as the region where the particles are injected;
  hvROI_p - a 1*4 vector recording the dynamic ROI on side view in last cycle, initialized as the region where the particles are injected;
  xyhv_cpre_mat - a 20*4 matrix recording position of the swarm in previous 20 cycles in pixels;
  Branch_p - a n*7 matrix recording the branch map from last cycle, n dynamically changes during exploration;
  PosSwm_p - a 1*4 vector showing position of the swarm in last cycle, [x,y,z,type] in millimeter;
  Region - region of exploration, [x_upperleft, y_upperleft, x_lowerright, y_lowerright] in pixels;
  Img_total_arr_p - exploration map from last cycle, resolution of 1400*900;
  ReconsData_p - a n*6 matrix recording the data for reconstruction from last cycle;
  pix_inc_p - a 1*4 vector recording the incremental area in pixels of the swarm in last cycle;
  num_frame - current number of image frames;
  frame_ct - counter for pause;
  mf_dir_p - yaw angle of the rotating magnetic field in last cycle. 

Output:
  PosSwm - a 1*4 vector showing position of the swarm, [x,y,z,type] in millimeter;
  xyhv_co_mat - a 20*4 matrix recording positions of the swarm till current frame;
  xyROI, hvROI - dynamic ROI on top view and side view, respectively;
  Branch - a n*7 matrix recording the updated branch map;
  Img_total_arr - a 1400*900 matrix recording the updated exploration map;
  ReconsData - a n*6 matrix recording the updated data for reconstruction;
  all_end - marker to complete exploration;
  pix_inc - a 1*4 vector recording the incremental area in pixels of the swarm;
  frame_ct - updated counter for pause;
  mf_dir_p - updataed yaw angle of the rotating magnetic field; 
  dist_hvroi - moving direction of the ROI on side view.

Instructions:
1. Write a loop in the main program of magnetic actuation system, and the loop should be parallel to the loops that controlling the magnetic actuation process; 
2. Initialize the vaiables according to the list above;
3. Acquire image frames from top-view camera and side-view camera, and both images are in the resolution of 1400*900;
4. Initialize the background side-view image as the first side-view image acquired;
5. Set the loop to execute every N frames, where N can be adjusted according to the demanded rate for image processing and decision-making;
6. Integrate the MATLAB function into the program (in LabVIEW, the "function" line in the beginning and the "end" line at the end of the function should be removed when using "MATLAB script");
7. Extract the ROI image with a resolution of 128*128 on top view using the vector that is about to serve as the input parameter "xyROI_p", convert the ROI image into greyscale and its resolution into 256*256, binarize the image, and assign the result to variable "imgroi2_top";
8. Extract the ROI image with a resolution of 128*400 on side view using the vector that is about to serve as the input parameter "hvROI_p", convert the ROI image into greyscale, and assign the result to variable "img_grey_side";
9. Extract the ROI image with a resolution of 128*400 on background side-view image using the vector that is about to serve as the input parameter "hvROI_p", convert the ROI image into greyscale, and assign the result to variable "img_grey_side_ini";
10. Set shift registers to store the values of all the output variables, so that their values can be used in the next cycle;
11. Update the background side-view image every time the ROI on side view changes moving directions according to the value of "dist_hvroi";
12. Convey the values of "PosSwm" and "mf_dir" to the loops controlling the magnetic actuation system to actuate the swarm.

Conducting active exploration demands the cooperation with the magnetic actuation platform, microrobot swarm and feebacks from imaging devices. Hardware designs of the magnetic actuation systems please refer to:
 - Du, X., Yang, L., Yu, J., Chan, K.F., Chiu, P.W.Y., Zhang, L.: Robomag: A magnetic actuation system based on mobile electromagnetic coils with tunable working space. In: 2020 5th International Conference on Advanced Robotics and Mechatronics (ICARM), pp. 125–131 (2020). IEEE
 - Du, X., Zhang, M., Yu, J., Yang, L., Chiu, P.W.Y., Zhang, L.: Design and real-time optimization for a magnetic actuation system with enhanced flexibility. IEEE/ASME Transactions on Mechatronics 26(3), 1524–1535 (2021)
 - Du, X., Yu, J.: Image-integrated magnetic actuation systems for localization and remote actuation of medical miniature robots: A survey. IEEE Transactions on Robotics 39(4), 2549–2568 (2023)
The principle and logic of the algorithm please refer to the Main Text, Methods and Supplementary Information in the manuscript. Demonstration of this function please refer to "Demo for Active Exploration\Demo from sample video.vi".
