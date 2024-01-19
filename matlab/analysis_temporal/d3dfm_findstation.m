function stationindex = d3dfm_findstation(hisfile,stationname,stationtype)
    
%read stations and create array

if any(strcmp(stationtype,{[],'','station'}))
    stations_raw=ncread(hisfile,'station_name');
    for n=1:size(stations_raw,2)
        stationArray{n,1}=strcat(stations_raw(:,n)');
    end
    %find index of to-be-analyzed station
    stationindex =find(strcmp(stationArray,stationname));
    
elseif any(strcmp(stationtype,{'cross section','cross_section','cross-section','xsection'}))
    cross_sections_raw=ncread(hisfile,'cross_section_name');
    for n=1:size(cross_sections_raw,2)
        crossSectionArray{n,1}=strcat(cross_sections_raw(:,n)');
    end
    %find index of to-be-analyzed station
    stationindex =find(strcmp(crossSectionArray,stationname));
    
else
    error('invalid type of station, choose between "station" or "cross-section"')
end



    
    
    