Copyright © 2023 DU Xingzhou.
Licence: Apache license 2.0.

Source code for real-time mapping, decision making and data recording in active exporation process. The code is in the format of MATLAB function and can be integrated to main program (LabView project in this work) when demanded. Conducting active exploration demands the cooperation with the magnetic actuation platform, microrobot swarm and feebacks from imaging devices. Please refer to our paper for further information. Demonstration of this function please refer to "Demo for Active Exploration\Demo from sample video.vi".

Input:
  imgroi2_top - binarized top-view image with a resolution of 1400*900;
  img_grey_side - greyscale side-view image with a resolution of 1400*900;
  img_grey_side_ini - greyscale background side-view image with a resolution of 1400*900;
  len_x, len_y - consts showing the dimension (in millimeter) of top view in x and y directions, respectively;
  ROI_len - length of the dynamic ROI on top view in pixels;
  ROI_hei - height of the dynamic ROI on side view in pixels;
  xyROI_p, hvROI_p - dynamic ROI on top view and side view in last cycle;
  xyhv_cpre_mat - position of the swarm in previous cycles;
  Branch_p - branch map from last cycle;
  PosSwm_p - position of the swarm in last cycle;
  Region - region of exploration;
  Img_total_arr_p - exploration map from last cycle;
  ReconsData_p - data for reconstruction from last cycle;
  pix_inc_p - incremental region of the swarm in last cycle;
  num_frame - current number of image frames;
  frame_ct - counter for pause;
  mf_dir_p - yaw angle of the rotating magnetic field in last cycle. 

Output:
  PosSwm - position of the swarm;
  xyhv_co_mat - positions of the swarm till current frame;
  xyROI, hvROI - dynamic ROI on top view and side view, respectively;
  Branch - updated branch map;
  Img_total_arr - updated exploration map;
  ReconsData - updated data for reconstruction;
  all_end - marker to complete exploration;
  pix_inc - incremental region of the swarm;
  frame_ct - updated counter for pause;
  mf_dir_p - updataed yaw angle of the rotating magnetic field; 
  dist_hvroi - moving direction of the ROI on side view.
