% mexPDallfaces.m Help file for mexPDall MEX-file
%  mexPDallfaces.cc - Uses the c++ library voro++ to generate 3D Laguerre diagrams
%  in cuboidal domains, with and without periodicity.
%
%  The standard calling syntax is
%
%     [v,t,xc,vfn]=mexPDallfaces(box,x,w,periodic)
%
%   where box is a 1x3 array [l1 l2 l3] specifying the container size
%         x   is an Nx3 array specifying the locations of the generators
%         w   is an Nx1 array specifying the weights
%         periodic is a boolean variable, if 'true' the computed diagram is
%         periodic if 'false' the computing diagram is not periodic
%
%   The output is
%         v     an Nx1 array of cell volumes
%         t     an Nx1 array of cell transport costs 
%         xc    an Nx3 array of cell centroids
%         vfn   an Nx3 cell array of cell vertex coordinates
%               vfn{i,1} contains an Nv x 3 array of vertex coordinates for the ith cell (Nv number of vertices)
%               vfn{i,2} contains an Nf x 3 cell array of vertex indices, face areas and outward face normals for each face of the ith cell (Nf number of faces)
%               vfn{i,3} contains an Nf x 1 array of neighbour indices for the ith cell
%   
%   If periodic is ommitted the assumption is that the periodic is 'false'
%   If w is ommitted then the assumption is that the weights are all zero
%   If box is ommitted then the assumption is that box=[1 1 1] so the
%   domain is the unit cube
%
%   It is therefore possible to call mexPDall with the following arguments
%
%     [v,t,xc,vfn]=mexPDallfaces(x);
%
%     [v,t,xc,vfn]=mexPDallfaces(box,x);
%     [v,t,xc,vfn]=mexPDallfaces(x,w);
%     [v,t,xc,vfn]=mexPDallfaces(x,periodic);
%
%     [v,t,xc,vfn]=mexPDallfaces(box,x,w);
%     [v,t,xc,vfn]=mexPDallfaces(box,x,periodic);
%     [v,t,xc,vfn]=mexPDallfaces(x,w,periodic);
