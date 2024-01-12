function d3dfm_mapdata=ncread_d3dfm_vector(mapfile,parameter_delft3dfm_nosuffix)
%   Description: 
%       Read x and y components of vector parameter, calculates magnitude and direction (degrees clockwise from North)
%       Function is used within "d3dfm_mapplot.m" 
%
%   Author: 
%       Martijn Bregman (created 7/12/2022)
%
%   Input:
%       mapfile: path to Delft3D FM netcdf map output file
%       parameter_delft3dfm_nosuffix: name of vector variable (as stored in netcdf output file), without the trailing "x" or "y"
%           e.g., when defining [parameter_delft3dfm_nosuffix='mesh2d_uc'], script will read 'mesh2d_ucx' and 'mesh2d_ucy'

% define parameter names in output netcdf file
parameter_delft3dfm_x=[parameter_delft3dfm_nosuffix 'x'];
parameter_delft3dfm_y=[parameter_delft3dfm_nosuffix 'y'];

% read x and y components
d3dfm_mapdata.x=ncread(mapfile,parameter_delft3dfm_x);
d3dfm_mapdata.y=ncread(mapfile,parameter_delft3dfm_y);

% compile magnitude and direction
d3dfm_mapdata.mag=sqrt(d3dfm_mapdata.x.^2+d3dfm_mapdata.y.^2);
d3dfm_mapdata.dir_degN=mod(90-atan2d(d3dfm_mapdata.y,d3dfm_mapdata.x),360) %y,x


