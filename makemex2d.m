
% makemex2d.m
%
% AUTHOR DB SR
% This script compilex the .cc files to build the mex-files mexPD_2d, mexPDall_2d

% Relative path to voro++ 2d source files
vorodir='voro_src/src_2d/';

% Source files
vorofiles={'cell_2d.cc' 'common.cc' 'container_2d.cc' 'v_base_2d.cc' 'v_compute_2d.cc' 'c_loops_2d.cc' 'wall_2d.cc' 'cell_nc_2d.cc' 'ctr_boundary_2d.cc' 'ctr_quad_2d.cc' 'quad_march.cc'};

eval([sprintf('mex -I%s ',vorodir) 'mexPD_2d.cc ' strjoin(strcat(vorodir,vorofiles))]);
eval([sprintf('mex -I%s ',vorodir) 'mexPDall_2d.cc ' strjoin(strcat(vorodir,vorofiles))]);
