function FMg = d3dfm_readmeshgeometry(datafile)
%   Description: 
%       Outputs a strurcture of model mesh information from DelftFM data file. 
%       Read coordinates for mesh nodes, mesh faces, edge nodes, and mesh relationships:
%          vertex nodes (counterclockwise), Neighboring faces of mesh edges
%       Performs a triangulaion of the face nodes, but does not output this currently
%
%       Function is used within "d3dfm_mapplot.m" 
%
%   Author: 
%       Chris Esposito
%
%   Input:
%        datafile: path to Delft3D FM netcdf map output file



%% nodes
% Coordinates of mesh nodes (including z)
% Dimensions: mesh2d_nNodes
% (number of nodes)
disp('nodes...')
node_x = ncread(datafile,'mesh2d_node_x');
node_y = ncread(datafile,'mesh2d_node_y');
node_z = ncread(datafile,'mesh2d_node_z');


%% faces
%Coordinates of mesh face
% Dimensions: mesh2d_nFaces
% 'Characteristic coordinate of mesh face'
disp('faces...')
face_x = ncread(datafile,'mesh2d_face_x');
face_y = ncread(datafile,'mesh2d_face_y');


%% edges
% edge nodes
% Dimensions: 2,mesh2d_nEdges
% 'Start and end nodes of mesh edges'
disp('edges...')
edge_nodes = ncread(datafile,'mesh2d_edge_nodes');

% edge midpoints for mesh
% Dimensions: mesh2d_nEdges
edge_x = ncread(datafile,'mesh2d_edge_x');
edge_y = ncread(datafile,'mesh2d_edge_y');

%% attributes
area=ncread(datafile,'mesh2d_flowelem_ba');


%% relationships
% Vertex nodes of mesh faces (counterclockwise)
%  Dimensions: mesh2d_nMax_face_nodes,mesh2d_nFaces   
% (max number of nodes on any face in the mesh, number of faces)
disp('relationships...')
face_nodes = ncread(datafile,'mesh2d_face_nodes');

% Dimensions: 2,mesh2d_nEdges
% 'Neighboring faces of mesh edges'
%%%%%%% I DO NOT KNOW WHICH FACE COMES FIRST %%%%%%%%
edge_faces = ncread(datafile,'mesh2d_edge_faces');      %why is this coming into matlab as a double? but in ncdisp it is listed as int32

%T is a triangulation, with larger shapes split into triangles
%Ti is the corresponding face indices for each triangle
%[T,Ti]=getTri(face_nodes);


%% everything in a struct
disp('all in a structure...')
FMg.node_x=node_x;
FMg.node_y=node_y;
FMg.node_z=node_z;
FMg.face_x=face_x;
FMg.face_y=face_y;
FMg.area=area;
FMg.edge_nodes=edge_nodes;
FMg.edge_x=edge_x;
FMg.edge_y=edge_y;
FMg.face_nodes=face_nodes;
FMg.edge_faces=edge_faces;