function d3dfm_mapplot(sim_outputpath,parameter_mag_modelname,varargin) 
%   Description: 
%       Plot maps of Delft3D FM output with optional averaging/differencing functionality, and vector overlay
%
%   Author: 
%       Martijn Bregman (created 1/5/2023)
%
%   Syntax:
%       d3dfm_mapplot(sim_outputpath,<keyword,value>)
%
%   Input:
%       sim_outputpath: path to simulation output folder
%       parameter_mag_modelname: name of variable (as stored in netcdf output file) 
%
%   Keyword-value pairs: 
%       see list w/examples below
%
%   Dependencies:
%       d3dfm_readmeshgeometry
%       d3dfm_readtime
%       d3dfm_plotmesh
%       d3dfm_processmapdata
%       ncread_d3dfm_vector
%



%% List of optional arguments (defined with default settings)
% basic processing options
processing.interval={'instantaneous'};  %{'instantaneous','daily','weekly','monthly','entiresimulation'}
processing.intervalProcessing='average'; % 'atBegin' (begin of interval), 'atEnd' (end of interval) 'average' (average per interval), 'change' (change), 'change_sinceTStart' (change since tstart)
processing.tstart=[]; % e.g. datenum(2020,01,01)
processing.tend=[]; 

% plotting options 
plotting.ldbfile=''; %e.g. 'sentinel2_LULC_southeastLA_crop.mat'
plotting.coordinatesystem=''; % e.g. 'UTM 15N' 
plotting.geographicareas_names={'default'};  % e.g. {'modeldomain','areaofinterest'} %use same value if defined once only, use array values is defined multiple times
plotting.easting_bounds={'auto'}; % e.g. {[750000 850000],[785000 810000]}
plotting.northing_bounds={'auto'}; % e.g. [3200000 3300000],[[3240000 3265000]]}
plotting.colormap_name='jet'; % e.g. 'turbo' 
plotting.colormap_limits={'auto'}; %e.g. [-1 3]
plotting.northarrow_filename=''; % e.g. 'WI_NorthArrow.png'
plotting.parameter_mag_printname=parameter_mag_modelname; % e.g. 'Water level (m NAVD88)'
plotting.parameter_mag_colorbarLabel='';


% Vector overlay (optional)
vectorplot.parameter_modelname=''
vectorplot.parameter_vec_printname=''; % e.g. 'Velocity vectors'
vectorplot.vector_count=''; % e.g. 10 (in each direction)
vectorplot.vector_scaling=''; %e.g. 0.5
vectorplot.scalevector_equivmagnitude=''; 
vectorplot.unit=''; %m/s'
vectorplot.inpolygon=''; %only display vectors within polygon

% output options
output.simname='';
output.destination_folder=pwd;
output.savefig=0;

%% parse optional inputs
for setting = 1:2:length(varargin)
    switch varargin{setting}

        case 'processing.interval'
            processing.interval = varargin{setting + 1};
        case 'processing.intervalProcessing'
            processing.intervalProcessing = varargin{setting + 1};
        case 'processing.tstart'
            processing.tstart = varargin{setting + 1};
        case 'processing.tend'
            processing.tend = varargin{setting + 1};  
            
        case 'plotting.ldbfile'
            plotting.ldbfile = varargin{setting + 1};  
        case 'plotting.coordinatesystem'
            plotting.coordinatesystem = varargin{setting + 1};  
        case 'plotting.geographicareas_names'
            plotting.geographicareas_names = varargin{setting + 1};  
        case 'plotting.easting_bounds'
            plotting.easting_bounds = varargin{setting + 1}; 
        case 'plotting.northing_bounds'
            plotting.northing_bounds = varargin{setting + 1};  
        case 'plotting.colormap_limits'
            plotting.colormap_limits = varargin{setting + 1}; 
        case 'plotting.colormap_name'
            plotting.colormap_name = varargin{setting + 1}; 
        case 'plotting.northarrow_filename'
            plotting.northarrow_filename = varargin{setting + 1};             
        case 'plotting.parameter_mag_printname' 
            plotting.parameter_mag_printname = varargin{setting + 1};             
        case 'plotting.parameter_mag_colorbarLabel' 
            plotting.parameter_mag_colorbarLabel = varargin{setting + 1};             
            
        case 'vectorplot.parameter_modelname'
            vectorplot.parameter_modelname = varargin{setting + 1}; 
        case 'vectorplot.parameter_vec_printname'
            vectorplot.parameter_vec_printname = varargin{setting + 1};  
        case 'vectorplot.vector_count' % #7 integer (in each direction)
            vectorplot.vector_count = varargin{setting + 1}; 
        case 'vectorplot.vector_scaling' % #8
            vectorplot.vector_scaling = varargin{setting + 1};            
        case 'vectorplot.scalevector_equivmagnitude'
            vectorplot.scalevector_equivmagnitude = varargin{setting + 1}; 
        case 'vectorplot.unit'
            vectorplot.unit = varargin{setting + 1}; 
        case 'vectorplot.inpolygon'
            vectorplot.inpolygon = varargin{setting + 1}; 
            
        case 'output.simname'
            output.simname = varargin{setting + 1};
        case 'output.destination_folder'
            output.destination_folder = varargin{setting + 1};
        case 'output.savefig'
            output.savefig = varargin{setting + 1};
        otherwise
            [];
    end
end

%% create output directory
if ~exist(output.destination_folder,'dir'); mkdir(output.destination_folder); end;
if ~exist([output.destination_folder filesep 'fig_files'],'dir') && output.savefig==1; mkdir([output.destination_folder filesep 'fig_files']); end;
% 
%% read data
if size(sim_outputpath,1)==1 %single simulation
    ncMap_dir=dir([sim_outputpath filesep '*map.nc'])
    for iMap=1:length(ncMap_dir)
        disp(['reading from map file ' num2str(iMap) '/' num2str(length(ncMap_dir))])
        tic
        FMgeometry{iMap,1}=d3dfm_readmeshgeometry([ncMap_dir(iMap,1).folder filesep ncMap_dir(iMap,1).name]); %organize data into mappable structure
        % read magnitude data
        d3dfm_mapdata_mag{iMap,1} = ncread([ncMap_dir(iMap,1).folder filesep ncMap_dir(iMap,1).name],parameter_mag_modelname);  %Data to show in faces
        if ndims(d3dfm_mapdata_mag{iMap,1})==3 %sum for all sediment fractions
            d3dfm_mapdata_mag{iMap,1}=squeeze(sum(d3dfm_mapdata_mag{iMap,1},2)); 
        end
        
        % read vector data if defined
        if ~isempty(vectorplot.parameter_modelname)
            %WIP: wrap next couple of lines in single function
            d3dfm_mapdata_vec{iMap,1} = ncread_d3dfm_vector([ncMap_dir(iMap,1).folder filesep ncMap_dir(iMap,1).name],vectorplot.parameter_modelname);  %Data to show in faces
            d3dfm_mapdata_vec_mag{iMap,1} = d3dfm_mapdata_vec{iMap,1}.mag;
            d3dfm_mapdata_vec_xcomp{iMap,1} = d3dfm_mapdata_vec{iMap,1}.x;
            d3dfm_mapdata_vec_ycomp{iMap,1} = d3dfm_mapdata_vec{iMap,1}.y;
            d3dfm_mapdata_vec_dir_degN{iMap,1} = d3dfm_mapdata_vec{iMap,1}.dir_degN;
        end
        toc
    end

% read two simulations and calculate difference    
elseif isempty(vectorplot.parameter_modelname) && size(sim_outputpath,1)==2 && size(sim_outputpath,2)==1 %difference between simulations
    ncMap_dir=dir([sim_outputpath{1,1} filesep '*map.nc'])
    ncMap_dir_ref=dir([sim_outputpath{2,1} filesep '*map.nc'])
    
    for iMap=1:length(ncMap_dir)
        disp(['reading from map file ' num2str(iMap) '/' num2str(length(ncMap_dir))])
        tic
        FMgeometry{iMap,1}=d3dfm_readmeshgeometry([ncMap_dir(iMap,1).folder filesep ncMap_dir(iMap,1).name]); %organize data into mappable structure
        % read magnitude data
        d3dfm_mapdata{iMap,1} =     ncread([ncMap_dir(iMap,1).folder filesep ncMap_dir(iMap,1).name],parameter_mag_modelname);  %Data to show in faces
        d3dfm_mapdata_ref{iMap,1} = ncread([ncMap_dir_ref(iMap,1).folder filesep ncMap_dir_ref(iMap,1).name],parameter_mag_modelname);  %Data to show in faces
        d3dfm_mapdata_mag{iMap,1} = d3dfm_mapdata{iMap,1}-d3dfm_mapdata_ref{iMap,1};
        
        if ndims(d3dfm_mapdata_mag{iMap,1})==3 %sum for all sediment fractions
            d3dfm_mapdata_mag{iMap,1}=squeeze(sum(d3dfm_mapdata_mag{iMap,1},2)); 
        end
 
        toc
    end
    
% elseif size(sim_outputpath,1)==2 && size(sim_outputpath,2)==2 %difference between change in two simulations (e.g. compare FWA vs FWOA bed level change, where initial topobathy was different)
%     ncMap_dir=dir([sim_outputpath{1,1} filesep '*map.nc'])
%     ncMap_dir_t0=dir([sim_outputpath{1,2} filesep '*map.nc'])
%     ncMap_dir_ref=dir([sim_outputpath{2,1} filesep '*map.nc'])
%     ncMap_dir_ref_t0=dir([sim_outputpath{2,2} filesep '*map.nc'])
%     %remove merged map if available
% %     ncMap_dir = ncMap_dir(~ismember({ncMap_dir.name},{'FlowFM_merged_map.nc'}));
% %     ncMap_dir_ref = ncMap_dir_ref(~ismember({ncMap_dir_ref.name},{'FlowFM_merged_map.nc'}));
% 
%     for iMap=1:length(ncMap_dir)
%         disp(['reading from map file ' num2str(iMap) '/' num2str(length(ncMap_dir))])
%         tic
%         FMgeometry{iMap,1}=d3dfm_readmeshgeometry([ncMap_dir(iMap,1).folder filesep ncMap_dir(iMap,1).name]); %organize data into mappable structure
%         d3dfm_mapdata{iMap,1} = ncread([ncMap_dir(iMap,1).folder filesep ncMap_dir(iMap,1).name],parameter_mag_modelname);  %Data to show in faces
%         d3dfm_mapdata_t0{iMap,1} = ncread([ncMap_dir_t0(iMap,1).folder filesep ncMap_dir_t0(iMap,1).name],parameter_mag_modelname);  %Data to show in faces
%         d3dfm_mapdata_change{iMap,1}=d3dfm_mapdata{iMap,1}-d3dfm_mapdata_t0{iMap,1};
%         
%         d3dfm_mapdata_ref{iMap,1} = ncread([ncMap_dir_ref(iMap,1).folder filesep ncMap_dir_ref(iMap,1).name],parameter_mag_modelname);  %Data to show in faces
%         d3dfm_mapdata_ref_t0{iMap,1} = ncread([ncMap_dir_ref_t0(iMap,1).folder filesep ncMap_dir_ref_t0(iMap,1).name],parameter_mag_modelname);  %Data to show in faces
%         d3dfm_mapdata_ref_change{iMap,1}=d3dfm_mapdata_ref{iMap,1}-d3dfm_mapdata_ref_t0{iMap,1};
%         
%         d3dfm_mapdata_mag{iMap,1}=d3dfm_mapdata_change{iMap,1}-d3dfm_mapdata_ref_change{iMap,1};
%         toc
%     end
else 
    disp('script cannot handle more than two simulations')
end
% 
% 
%% read time
iMap=1;
model_time=d3dfm_readtime([ncMap_dir(iMap,1).folder filesep ncMap_dir(iMap,1).name]);
% find indices of model data closest to start and end of period under consideration
if ~isempty(processing.tstart) 
    [~,plotting.tstart_index]=min(abs(model_time-processing.tstart));
else
    plotting.tstart_index=1; %start of simulation 
end
if ~isempty(processing.tend) 
    [~,plotting.tend_index]=min(abs(model_time-processing.tend));
else
    plotting.tend_index=length(model_time); %start of simulation 
end

% read land boundary file
if ~isempty(plotting.ldbfile); 
    if strfind(plotting.ldbfile,'.mat')
        load(plotting.ldbfile);
    elseif strfind(plotting.ldbfile,'.ldb')
        plotting.ldbdata=readldb(plotting.ldbfile); 
    end
end
%% average and plot for each of the selected averaging options

if ~iscell(processing.interval)
    error('"processing.interval" has be to be defined as cell array')
end
for i_interval=1:length(processing.interval)
    interval=processing.interval{i_interval};  % define averaging option
    % average data    
    [d3dfm_mapdata_mag_processed,processing_namearray]=d3dfm_processmapdata(d3dfm_mapdata_mag,model_time,interval,processing.intervalProcessing,plotting.tstart_index,plotting.tend_index);  

    if ~isempty(vectorplot.parameter_modelname)
        %WIP: wrap vector averaging and concatenating into function        
        %average
        [d3dfm_mapdata_vec_mag_processed,~]=d3dfm_processmapdata(d3dfm_mapdata_vec_mag,model_time,interval,processing.intervalProcessing,plotting.tstart_index,plotting.tend_index);  
        [d3dfm_mapdata_vec_xcomp_processed,~]=d3dfm_processmapdata(d3dfm_mapdata_vec_xcomp,model_time,interval,processing.intervalProcessing,plotting.tstart_index,plotting.tend_index);      
        [d3dfm_mapdata_vec_ycomp_processed,~]=d3dfm_processmapdata(d3dfm_mapdata_vec_ycomp,model_time,interval,processing.intervalProcessing,plotting.tstart_index,plotting.tend_index);  
        
        %compile concatenated arrays
        d3dfm_mapdata_face_x=[];
        d3dfm_mapdata_face_y=[];
        for iMap=1:length(ncMap_dir)
            d3dfm_mapdata_vec_dir_degN_processed{iMap,1}=mod(90-atan2d(d3dfm_mapdata_vec_ycomp_processed{iMap,1},d3dfm_mapdata_vec_xcomp_processed{iMap,1}),360); %y,x
            d3dfm_mapdata_face_x=[d3dfm_mapdata_face_x; FMgeometry{iMap,1}.face_x];
            d3dfm_mapdata_face_y=[d3dfm_mapdata_face_y; FMgeometry{iMap,1}.face_y];
        end
        d3dfm_mapdata_vec_xcomp_processed_cat=cat(1,d3dfm_mapdata_vec_xcomp_processed{:});
        d3dfm_mapdata_vec_ycomp_processed_cat=cat(1,d3dfm_mapdata_vec_ycomp_processed{:});
        d3dfm_mapdata_vec_mag_processed_cat=cat(1,d3dfm_mapdata_vec_mag_processed{:});
        d3dfm_mapdata_vec_dir_degN_processed_cat=cat(1,d3dfm_mapdata_vec_dir_degN_processed{:});
        
    end

    %% plotting  
    if ~iscell(plotting.geographicareas_names)
    error('"plotting.geographicareas_names" has be to be defined as cell array')
    end
    for itimestep=1:length(processing_namearray)
        plotting_figure=figure; hold on;
        set(plotting_figure,'Position',[0 0 900 600]);
        % plot output
        for iMap=1:length(ncMap_dir)
            d3dfm_plotmesh(FMgeometry{iMap,1},'facedata',d3dfm_mapdata_mag_processed{iMap,1}(:,itimestep));
        end
        % plot land boundary
        if ~isempty(plotting.ldbfile)
            plot(plotting.ldbdata.x,plotting.ldbdata.y,'color',[0.500    0.500    0.500],'linewidth',0.5)
        end
        box on; grid on; axis equal;

        xlabel(['Easting (' plotting.coordinatesystem ' km)']);
        ylabel(['Northing (' plotting.coordinatesystem ' km)']);
        
        for iplot=1:length(plotting.geographicareas_names)
           
            plotting.charthandle=gca;
            plotting.charthandle.LineWidth = 2; % increase box line thickness
            % delete north arrow if existent
            axesHandlesToChildObjects = findobj(plotting.charthandle, 'Type', 'image');
            if ~isempty(axesHandlesToChildObjects)
                delete(axesHandlesToChildObjects);
            end
            % define properties based on user input
            plotname=plotting.geographicareas_names{iplot};
            
            xlim(plotting.easting_bounds{iplot}); ylim(plotting.northing_bounds{iplot});
            caxis(plotting.colormap_limits{iplot});
            colormap(plotting.colormap_name); colorbar;
            ylabel(colorbar,plotting.parameter_mag_colorbarLabel,'FontSize',11)
            % scale ticks to kilometers
            plot_ticks_x = num2cell(round(get(plotting.charthandle,'xtick')/1000));
            plot_ticks_y = num2cell(round(get(plotting.charthandle,'ytick')/1000));
            set(plotting.charthandle,'xticklabel',plot_ticks_x,'yticklabel',plot_ticks_y)
            
            if ~isempty(vectorplot.parameter_modelname) && ~isempty(vectorplot.vector_count)
                % set up grid for interpolation
                charthandle=gca;
                plot_xlim = get(charthandle,'xlim');
                plot_ylim = get(charthandle,'ylim');
                %vectorplot_dx=600; %distance between grid points
                vectorplot_dx=mean([diff(plot_ylim)/vectorplot.vector_count,diff(plot_xlim)/vectorplot.vector_count]); %dynamically determine distance between grid points
                vectorplot_xpoints = plot_xlim(1)+0.5*vectorplot_dx:vectorplot_dx:plot_xlim(2)-0.5*vectorplot_dx; %grid x points, start and end from axes lims
                vectorplot_ypoints = plot_ylim(1)+0.5*vectorplot_dx:vectorplot_dx:plot_ylim(2)-0.5*vectorplot_dx; %grid y points
                [vectorplot_xmesh,vectorplot_ymesh] = meshgrid(vectorplot_xpoints,vectorplot_ypoints);
                
                % interpolate to grid
                interpolant_xcomp = scatteredInterpolant(d3dfm_mapdata_face_x, d3dfm_mapdata_face_y,d3dfm_mapdata_vec_xcomp_processed_cat(:,itimestep));
                interpolant_ycomp = scatteredInterpolant(d3dfm_mapdata_face_x, d3dfm_mapdata_face_y,d3dfm_mapdata_vec_ycomp_processed_cat(:,itimestep));
                vectorplot_x_int = interpolant_xcomp(vectorplot_xmesh,vectorplot_ymesh);
                vectorplot_y_int = interpolant_ycomp(vectorplot_xmesh,vectorplot_ymesh);
                
                % remove interpolated data outside of grid perimeter
                if ~isempty(vectorplot.inpolygon)
                    vectorplot_perimeter=readldb(vectorplot.inpolygon);
                    vectorplot_insideperimeter=inpolygon(vectorplot_xmesh,vectorplot_ymesh,vectorplot_perimeter.x,vectorplot_perimeter.y);
                    vectorplot_x_int(~vectorplot_insideperimeter)=NaN;
                    vectorplot_y_int(~vectorplot_insideperimeter)=NaN;
                end


                
                % ADD SCALE ARROW
                if ~isempty(vectorplot.scalevector_equivmagnitude)
                    % define number of indices that need to be ommitted to make room for scale arrow
%                     vectorplot_scalearrow_verticalmargin=ceil(0.035*vectorplot.vector_count); % omit top 12%
                    vectorplot_scalearrow_verticalmargin=ceil(0.10*vectorplot.vector_count); % omit top 12%
                    vectorplot_scalearrow_horizontalmargin=ceil(0.35*vectorplot.vector_count); % omit left 35%
                    
                    % replace values with 0 in a given margin from the top left of the plot (bottom left in table)
                    vectorplot_xindices_zeroed=...
                        [size(vectorplot_x_int,1)+1-vectorplot_scalearrow_verticalmargin]:... %some number of indices before the end
                        size(vectorplot_x_int,1); %up to the end
                    vectorplot_yindices_zeroed=1:vectorplot_scalearrow_horizontalmargin; %up to some number from the start
                    vectorplot_x_int(vectorplot_xindices_zeroed,vectorplot_yindices_zeroed)=0;
                    vectorplot_y_int(vectorplot_xindices_zeroed,vectorplot_yindices_zeroed)=0;
                    
                    % define magnitude of scale vector (located in top left of plot = bottom left in matrix)
                    vectorplot_x_int(size(vectorplot_y_int,1),1)=vectorplot.scalevector_equivmagnitude;
                    
                    % define position of annotion (small margin if there's only a single-index gap, or two indices below when there's a larger gap
                    if vectorplot_scalearrow_verticalmargin==1
                        yposition=vectorplot_ypoints(size(vectorplot_x_int,1))-0.2*vectorplot_dx;
                    else
                        yposition=vectorplot_ypoints(size(vectorplot_x_int,1)-2);
                    end
                    % add annotation to vector
                    vectorplot_arrowann=text(vectorplot_xpoints(1),yposition,...
                        ['Arrow scale = ' num2str(vectorplot.scalevector_equivmagnitude) ' ' vectorplot.unit],...
                        'BackgroundColor','w','EdgeColor','k');
                    %  delete([vectorplot_base,vectorplot_overlay,vectorplot_arrowann]) % for testing only
                end
                % plot vectors (thick white arrow with thin black overlay for contrast)
                vectorplot_base=quiver(charthandle,vectorplot_xmesh,vectorplot_ymesh,vectorplot_x_int, vectorplot_y_int,1.1*vectorplot.vector_scaling,'color','w','LineWidth',2);
                vectorplot_overlay=quiver(charthandle,vectorplot_xmesh,vectorplot_ymesh,vectorplot_x_int, vectorplot_y_int,vectorplot.vector_scaling,'color','k','LineWidth',0.5); %,'MarkerSize',1,'MaxHeadSize',10)

                % add title
                if strcmp(processing_namearray{itimestep,1},'singledt')
                    timestepname='';
                else
                    timestepname=[' | ' processing_namearray{itimestep,1}];
                end
                title({[plotting.parameter_mag_printname ' | ' vectorplot.parameter_vec_printname timestepname],''},'Interpreter','none');
            else
                if strcmp(processing_namearray{itimestep,1},'singledt'); 
                    timestepname='';
                else
                    timestepname=[' | ' processing_namearray{itimestep,1}];
                end
                title({[plotting.parameter_mag_printname timestepname],''},'Interpreter','none');

            end
            
            % add north arrow, position in upper right based on figure limits
            if ~isempty(plotting.northarrow_filename)
                [plotting_north_arrow,~,alpha] = imread(plotting.northarrow_filename);
                plotting_north_arrow = flipud(plotting_north_arrow);
                XLIM=get(gca, 'XLim');
                YLIM=get(gca, 'YLim');
                plotting_north_arrow_x_size = diff(XLIM)*0.08; % 8% of the width (roughly twice as tall as it is wide)
                plotting_north_arrow_y_size = diff(YLIM)*0.14; % 14% of the hieght
                plotting_north_arrow_x_buffer = diff(XLIM)*0.01; %buffer so it is not right on the edge
                plotting_north_arrow_y_buffer = diff(YLIM)*0.01;
                plotting_north_arrow_x_pos = [XLIM(2)-(plotting_north_arrow_x_size+plotting_north_arrow_x_buffer) XLIM(2)-plotting_north_arrow_x_buffer];
                plotting_north_arrow_y_pos = [YLIM(2)-(plotting_north_arrow_y_size+plotting_north_arrow_y_buffer) YLIM(2)-plotting_north_arrow_y_buffer];
                %display north arrow
                plotting_north_arrow_image=image(plotting.charthandle,'CData',plotting_north_arrow,'XData',plotting_north_arrow_x_pos,'YData',plotting_north_arrow_y_pos);
                if contains(plotting.northarrow_filename,'inverted')
                    temp = sum(plotting_north_arrow,3)/765; %make background transparent
                else
                    temp = 1 - sum(plotting_north_arrow,3)/765; %make background transparent
                end
                plotting_north_arrow_image.AlphaData = temp;
            end
            % save output
            plotting.figurename=[output.simname '__' parameter_mag_modelname '__' plotname '__' processing_namearray{itimestep,1}];
            print(plotting_figure,'-dpng', '-r150',[output.destination_folder filesep plotting.figurename '.png']);
            if output.savefig==1
                savefig(plotting_figure,[output.destination_folder filesep 'fig_files' filesep plotting.figurename '.fig'])
            end
            % delete vectors and arrow for next plot
            if exist('vectorplot_arrowann'); delete(vectorplot_arrowann); end;
            if exist('vectorplot_base'); delete(vectorplot_base); end;
            if exist('vectorplot_overlay'); delete(vectorplot_overlay); end;

        end
        close(plotting_figure)
    end

end
disp('d3dfm_mapplot finished')