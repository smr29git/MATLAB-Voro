
% makemex.m
%
% AUTHOR DB SR
% This script compiles the .cc files to build the three mex files mexPD, mexPDall, mexPDallfaces

% Relative path to voro++ sources files
vorodir='voro_src/src_3d/';

% Source files
vorofiles={'cell.cc' 'common.cc' 'container.cc' 'unitcell.cc' 'v_compute.cc' 'c_loops.cc' 'v_base.cc' 'wall.cc' 'pre_container.cc' 'container_prd.cc'};

% Build mexPD
eval(['mex mexPD.cc ' strjoin(strcat(vorodir,vorofiles))]);

% Build mexPDall
eval(['mex mexPDall.cc ' strjoin(strcat(vorodir,vorofiles))]);

% Build mexPDallfaces
eval(['mex mexPDallfaces.cc ' strjoin(strcat(vorodir,vorofiles))]);

