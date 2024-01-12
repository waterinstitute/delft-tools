function model_time = d3dfm_readtime(mapfile)
%   Description: 
%       Read time information from DelftFM map data file. 
%       "model_time" is returned in serial date numbers (e.g. similar to MATLAB's datenum)
%
%   Author: 
%       Martijn Bregman (created 6/8/2022)
%
%   Input:
%        mapfile: path to Delft3D FM netcdf map output file


model_timestr = strsplit(ncreadatt(mapfile,'time','units'),' ');% starting time of simulation
model_timeref = datenum([model_timestr{3} ' ' model_timestr{4}]); %starting time in datenumber
model_timeinsec = ncread(mapfile,'time'); % output time
model_time = model_timeref + model_timeinsec/3600/24;
