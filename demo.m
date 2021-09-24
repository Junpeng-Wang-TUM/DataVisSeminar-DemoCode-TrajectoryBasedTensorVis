%%Demonstration Code for Data Visualization Seminar Offered by the Chair of Computer Graphics and Visualization at TUM
%%Seminor Topic: Trajectory-based Tensor Field Visualization
%%Author: Junpeng Wang (junpeng.wang@tum.de)
%%Date: 2021-09-14
clear
clc

stressfileName = './data/cantilever2D_R500_iLoad5.vtk';
PSLsDensityCtrl = 20; %% control the density of the stress trajectories (positively-related)

TSV2D(stressfileName,PSLsDensityCtrl);
