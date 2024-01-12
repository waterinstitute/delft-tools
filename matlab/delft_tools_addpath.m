%   Description: script to add functions to MATLAB path
%
%   Author: Martijn Bregman (mbregman@thewaterinstitute.org)
%
%   Run script to add delft-tools to search path
%
%   Created: 1/10/2024
%

basepath = fileparts(mfilename('fullpath'));
filelist = dir(fullfile(basepath, '**\*.*'));  %get list of files and folders in any subfolder
filelist = filelist([filelist.isdir]);  %remove files from list, only keep folders

for ifile=1:length(filelist)
   addpath([filelist(ifile).folder filesep filelist(ifile).name]) 
end

disp('delft-tools successfully added to search path');