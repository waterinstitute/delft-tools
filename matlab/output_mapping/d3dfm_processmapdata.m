function [d3dfm_mapdata_processed,processing_namearray]=d3dfm_processmapdata(d3dfm_mapdata,model_time,interval,intervalProcessing,processing_t_start_ind,processing_t_end_ind)
%   Description: 
%       Time-average Delft3D FM map data. Function is used within "d3dfm_mapplot.m" 
%
%   Author: 
%       Martijn Bregman (created 6/8/2022)
%
%   Input:
%        d3dfm_mapdata: 2D array representing model map data
%        model_time: 1D array of 'datenum' values for model timestamps (generated with d3dfm_readtime)
%        interval: Timespan for data processing, choose from 'instantaneous', 'daily', 'weekly', 'monthly', or 'entiresimulation'
%        intervalProcessing: Approach for each interval - choose 'atBegin' (data at the interval's start), 'atEnd' (data at the interval's end), 'average' (average per interval), 'change' (change within the interval), 'change_sinceTStart' (change since processing_t_start_ind)
%        processing_t_start_ind: Start 'datenum' of the timeframe to be processed
%        processing_t_end_ind: End 'datenum' of the timeframe to be processed




%% time averaging
if strcmp(interval,'monthly')
    clear FMdata_averaged processing_namearray
    model_time_monthnumber=month(datetime(model_time,'ConvertFrom','datenum'));
    model_time_monthnumber_unique=model_time_monthnumber(diff([0;model_time_monthnumber])~=0);
    
    % define start and end month number, and array covering start month through end month
    processing_t_start_month=model_time_monthnumber(processing_t_start_ind);
    processing_t_end_month=model_time_monthnumber(processing_t_end_ind);
    if processing_t_end_month<processing_t_start_month
        processing_month_array=[processing_t_start_month:max(model_time_monthnumber) min(model_time_monthnumber):processing_t_end_month];
    else
        processing_month_array=[processing_t_start_month:processing_t_end_month];
    end
    
    % loop over months
    for imonth=1:length(processing_month_array)
        % find matching indices of model output for weeks under consideration
        processing_month_matching_indices{imonth,1}=find(ismember(model_time_monthnumber,processing_month_array(imonth)));
        % find calendar year associated with each week for naming purposes
        processing_month_matching_year(imonth,1)=year(datetime(datevec(model_time(processing_month_matching_indices{imonth,1}(1)))));
        % create name label from year and week numbers
        processing_namearray{imonth,1}=['monthly_' num2str(processing_month_matching_year(imonth,1)) '_' sprintf('%02.f', processing_month_array(imonth)) '_' intervalProcessing];
        
        if isequal(intervalProcessing,'atBegin')
            % pick first timestep
            for iMap=1:length(d3dfm_mapdata)
                d3dfm_mapdata_processed{iMap,1}(:,imonth)=d3dfm_mapdata{iMap,1}(:,processing_month_matching_indices{imonth,1}(1));
            end
        elseif isequal(intervalProcessing,'atEnd')
            % pick first timestep
            for iMap=1:length(d3dfm_mapdata)
                d3dfm_mapdata_processed{iMap,1}(:,imonth)=d3dfm_mapdata{iMap,1}(:,processing_month_matching_indices{imonth,1}(end));
            end            
        elseif isequal(intervalProcessing,'average')
            % average data
            for iMap=1:length(d3dfm_mapdata)
                d3dfm_mapdata_processed{iMap,1}(:,imonth)=mean(d3dfm_mapdata{iMap,1}(:,processing_month_matching_indices{imonth,1}),2);
            end
        elseif isequal(intervalProcessing,'change_sinceTStart')
            for iMap=1:length(d3dfm_mapdata)
                d3dfm_mapdata_processed{iMap,1}(:,imonth)=d3dfm_mapdata{iMap,1}(:,processing_month_matching_indices{imonth,1}(end))...
                    -d3dfm_mapdata{iMap,1}(:,processing_t_start_ind);                    
            end
        elseif isequal(intervalProcessing,'change')
            for iMap=1:length(d3dfm_mapdata)
                d3dfm_mapdata_processed{iMap,1}(:,imonth)=d3dfm_mapdata{iMap,1}(:,processing_month_matching_indices{imonth,1}(end))...
                    -d3dfm_mapdata{iMap,1}(:,max(processing_month_matching_indices{imonth,1}(1)-1,1));                    
            end                        
        else 
            disp('invalid mode selected for processing data')
        end
    end
    
elseif strcmp(interval,'weekly')
    clear FMdata_averaged processing_namearray
    % define week numbers of year
    %model_time_weekofyear=week(datetime(model_time,'ConvertFrom','datenum'),'weekofyear');
    % inconsistent numbering in default matlab function
    model_time_weekofyear=weeknumber(model_time)
    model_time_weekofyear_unique=model_time_weekofyear(diff([0;model_time_weekofyear])~=0);
    
    % define start and end week number, and array covering start week through end week
    processing_t_start_week=model_time_weekofyear(processing_t_start_ind);
    processing_t_end_week=model_time_weekofyear(processing_t_end_ind);
    if processing_t_end_week<processing_t_start_week
        processing_week_array=[processing_t_start_week:max(model_time_weekofyear) min(model_time_weekofyear):processing_t_end_week];
    else
        processing_week_array=[processing_t_start_week:processing_t_end_week];
    end
    
    % loop over weeks
    for iweek=1:length(processing_week_array)
        % find matching indices of model output for weeks under consideration
        processing_week_matching_indices{iweek,1}=find(ismember(model_time_weekofyear,processing_week_array(iweek)));
        % find calendar year associated with each week for naming purposes
        processing_week_matching_year(iweek,1)=year(datetime(datevec(model_time(processing_week_matching_indices{iweek,1}(1)))));
        % create name label from year and week numbers
        processing_namearray{iweek,1}=['weekly_' num2str(processing_week_matching_year(iweek,1)) '_' sprintf('%02.f', processing_week_array(iweek)) '_' intervalProcessing];

        if isequal(intervalProcessing,'atBegin')
            % pick first timestep
            for iMap=1:length(d3dfm_mapdata)
                d3dfm_mapdata_processed{iMap,1}(:,iweek)=d3dfm_mapdata{iMap,1}(:,processing_week_matching_indices{iweek,1}(1));
            end
        elseif isequal(intervalProcessing,'atEnd')
            % pick first timestep
            for iMap=1:length(d3dfm_mapdata)
                d3dfm_mapdata_processed{iMap,1}(:,iweek)=d3dfm_mapdata{iMap,1}(:,processing_week_matching_indices{iweek,1}(end));
            end            
        elseif isequal(intervalProcessing,'average')
        % average data
            for iMap=1:length(d3dfm_mapdata)
               d3dfm_mapdata_processed{iMap,1}(:,iweek)=mean(d3dfm_mapdata{iMap,1}(:,processing_week_matching_indices{iweek,1}),2);
            end
        elseif isequal(intervalProcessing,'change_sinceTStart')
            for iMap=1:length(d3dfm_mapdata)
                d3dfm_mapdata_processed{iMap,1}(:,iweek)=d3dfm_mapdata{iMap,1}(:,processing_week_matching_indices{iweek,1}(end))...
                    -d3dfm_mapdata{iMap,1}(:,processing_t_start_ind);                    
            end
        elseif isequal(intervalProcessing,'change')
            for iMap=1:length(d3dfm_mapdata)
                d3dfm_mapdata_processed{iMap,1}(:,iweek)=d3dfm_mapdata{iMap,1}(:,processing_week_matching_indices{iweek,1}(end))...
                    -d3dfm_mapdata{iMap,1}(:,max(processing_week_matching_indices{iweek,1}(1)-1,1));                    
            end              
        else 
            disp('invalid mode selected for processing data')
        end
    end
    
elseif strcmp(interval,'daily')
    clear FMdata_averaged processing_namearray
    % define day numbers of year
    model_time_dayofyearnumber=day(datetime(model_time,'ConvertFrom','datenum'),'dayofyear');
    
    % define start and end day number, and array covering start day through end day
    processing_t_start_day=model_time_dayofyearnumber(processing_t_start_ind);
    processing_t_end_day=model_time_dayofyearnumber(processing_t_end_ind);
    if processing_t_end_day<processing_t_start_day
        processing_day_array=[processing_t_start_day:max(model_time_dayofyearnumber) min(model_time_dayofyearnumber):processing_t_end_day];
    else
        processing_day_array=[processing_t_start_day:processing_t_end_day];
    end
    
    % loop over days
    for iday=1:length(processing_day_array)
        % find matching indices of model output for days under consideration
        processing_day_matching_indices{iday,1}=find(ismember(model_time_dayofyearnumber,processing_day_array(iday)));
        % obtain calendar year, month number, day of month number, associated with each day for naming purposes
        processing_day_matching_year(iday,1)=year(datetime(datevec(model_time(processing_day_matching_indices{iday,1}(1)))));
        processing_day_matching_month(iday,1)=month(datetime(datevec(model_time(processing_day_matching_indices{iday,1}(1)))));
        processing_day_matching_dayofmonth(iday,1)=day(datetime(datevec(model_time(processing_day_matching_indices{iday,1}(1)))));
        
        % create name label from year and week numbers
        processing_namearray{iday,1}=['daily_' num2str(processing_day_matching_year(iday,1)) '_' sprintf('%02.f', processing_day_matching_month(iday,1)) '_' sprintf('%02.f', processing_day_matching_dayofmonth(iday)) '_' intervalProcessing];

        if isequal(intervalProcessing,'atBegin')
        % pick first timestep
            for iMap=1:length(d3dfm_mapdata)
               d3dfm_mapdata_processed{iMap,1}(:,iday)=d3dfm_mapdata{iMap,1}(:,processing_day_matching_indices{iday,1}(1));
            end        
        elseif isequal(intervalProcessing,'atEnd')
        % pick first timestep
            for iMap=1:length(d3dfm_mapdata)
               d3dfm_mapdata_processed{iMap,1}(:,iday)=d3dfm_mapdata{iMap,1}(:,processing_day_matching_indices{iday,1}(end));
            end
        elseif isequal(intervalProcessing,'average')
        % average data
            for iMap=1:length(d3dfm_mapdata)
               d3dfm_mapdata_processed{iMap,1}(:,iday)=mean(d3dfm_mapdata{iMap,1}(:,processing_day_matching_indices{iday,1}),2);
            end
        elseif isequal(intervalProcessing,'change_sinceTStart')
            for iMap=1:length(d3dfm_mapdata)
                d3dfm_mapdata_processed{iMap,1}(:,iday)=d3dfm_mapdata{iMap,1}(:,processing_day_matching_indices{iday,1}(end))...
                    -d3dfm_mapdata{iMap,1}(:,processing_t_start_ind);
            end
        elseif isequal(intervalProcessing,'change')
            for iMap=1:length(d3dfm_mapdata)
                d3dfm_mapdata_processed{iMap,1}(:,iday)=d3dfm_mapdata{iMap,1}(:,processing_day_matching_indices{iday,1}(end))...
                    -d3dfm_mapdata{iMap,1}(:,max(processing_day_matching_indices{iday,1}(1)-1,1));                    
            end               
        else 
            disp('invalid mode selected for processing data')
        end       
    end
    
elseif strcmp(interval,'instantaneous')
    clear FMdata_averaged processing_namearray
    % loop over all output timestep, rearrange data without any averaging
    ctr=1;
    for model_itimestep=processing_t_start_ind:processing_t_end_ind
        
        if isequal(intervalProcessing,'average') || isequal(intervalProcessing,'atBegin') || isequal(intervalProcessing,'atEnd')
        % average data
            for iMap=1:length(d3dfm_mapdata)
                d3dfm_mapdata_processed{iMap,1}(:,ctr)=d3dfm_mapdata{iMap,1}(:,model_itimestep);
            end
            processing_namearray{ctr,1}=[datestr(model_time(model_itimestep),30)];
        elseif isequal(intervalProcessing,'change_sinceTStart')
            for iMap=1:length(d3dfm_mapdata)
                d3dfm_mapdata_processed{iMap,1}(:,ctr)=d3dfm_mapdata{iMap,1}(:,model_itimestep)-d3dfm_mapdata{iMap,1}(:,processing_t_start_ind);
            end
            processing_namearray{ctr,1}=[datestr(model_time(model_itimestep),30) '_' intervalProcessing];
        elseif isequal(intervalProcessing,'change')
            for iMap=1:length(d3dfm_mapdata)
                d3dfm_mapdata_processed{iMap,1}(:,ctr)=d3dfm_mapdata{iMap,1}(:,model_itimestep)-d3dfm_mapdata{iMap,1}(:,max(model_itimestep-1),1);
            end              
            processing_namearray{ctr,1}=[datestr(model_time(model_itimestep),30) '_' intervalProcessing];
        else 
            disp('invalid mode selected for processing data')
        end
  
        
        ctr=ctr+1;
    end
    
elseif strcmp(interval,'entiresimulation')
    clear FMdata_averaged processing_namearray

    processing_namearray{1,1}=['entiresimulation_' intervalProcessing];
    if isequal(intervalProcessing,'atBegin')
        % pick first timestep
        for iMap=1:length(d3dfm_mapdata)
            d3dfm_mapdata_processed{iMap,1}(:,1)=d3dfm_mapdata{iMap,1}(:,processing_t_start_ind);
        end
    elseif isequal(intervalProcessing,'atEnd')
        % pick first timestep
        for iMap=1:length(d3dfm_mapdata)
            d3dfm_mapdata_processed{iMap,1}(:,1)=d3dfm_mapdata{iMap,1}(:,processing_t_end_ind);
        end        
    elseif isequal(intervalProcessing,'average')
        % average data
        for iMap=1:length(d3dfm_mapdata)
            d3dfm_mapdata_processed{iMap,1}(:,1)=mean(d3dfm_mapdata{iMap,1}(:,1:processing_t_end_ind),2);
        end
    elseif isequal(intervalProcessing,'change_sinceTStart') || isequal(intervalProcessing,'change')
        for iMap=1:length(d3dfm_mapdata)
            d3dfm_mapdata_processed{iMap,1}(:,1)=d3dfm_mapdata{iMap,1}(:,processing_t_end_ind)-d3dfm_mapdata{iMap,1}(:,processing_t_start_ind);
        end
    else
        disp('invalid mode selected for processing data')
    end

    
elseif strcmp(interval,'singledt')
    clear FMdata_averaged processing_namearray
    for iMap=1:length(d3dfm_mapdata)
        d3dfm_mapdata_processed{iMap,1}=d3dfm_mapdata{iMap,1};
    end
    processing_namearray{1,1}='singledt';    
end