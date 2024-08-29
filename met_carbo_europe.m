%% met_carbo_europe.m
% check and process the met data from Dripsey/Wexford/Kerry for Carbo Europe
% for all years *not just 2005*.
%
% PL 1/7/2004
% PL modified 12.2005 to use new data structs.
% PL modified 21.08.2006. Bug in replacement of NaNs with 9999. Added more
% filters.
% PL 19.11.2007. Added PAR, Rswin zero error correction.
% 13.02.2009. Converted to function. Renamed to met_carbo_europe.m

% [1]   read in the met data
% [2]   select a row range ( = time coverage)
% [3]   manipluation:
%    change 0000 to 2400 if necessary
%    pressure -> kPa
% [4   swap columns (to get correct carbo europe order
% [5]   do checking
% [6]   write to output file ready for carbo europe. set all NaNs to -9999
%
% INPUTS
% met   met structure
% s     site parameters strcture
%
% OUTPUTS
%
function result=met_carbo_europe(met,s)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  general settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
startup; % path names

sitename = s.name; %
% define time range of interest
jd_start = 1;   % from AND? including
jd_end = 366;   % up to AND including
%ROOT_DIR = 'f:\Data\CARBOEUROPE\';
MISSING_VALUE=-9999;
time_format = '%10s\t%02d:%02d\t';  % this part of format string the same for all sites
header_file=fullfile(PATH_CARBOEUROPE,'met_header.txt');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  end of general settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% [1] read  data. set site specific parameters.
 switch upper(sitename{1})
     case 'IE-WEX'
         met.propwindspeed=met.two_d_winddir;
         met.propwinddir=met.two_d_windspeed;
 end
%
%         num_cols=29;
%         full_sitename='IE_DRI';
%         %INFILE = 'E:\Data\Donoughmore\30min_co2\raw_data\2003\T0_1_365_03fixed.dat';
%         %INFILE = 'E:\Data\Donoughmore\30min_co2\raw_data\2004\T0_1_366_04.dat';
%         %INFILE = 'f:\Data\Donoughmore\30min_co2\raw_data\2005\T0_1_365_05.dat';
%         % INFILE = 'f:\Data\Donoughmore\30min_co2\raw_data\2006\T0_1_365_06.dat';
%         %INFILE='F:\data\CARBOEUROPE\MQC_1_365_06.mat';
%         INFILE='f:\Data\Donoughmore\30min_co2\raw_data\2007\TO_1_274_07.csv';
%
%         % postprocess:
%         [f met s]=read_cork_data_new(INFILE,jd_start,jd_end);
%         %load(INFILE);
%         num_cols=29; % no of cols in output matrix
%         format_str = [time_format ...
%             '%6.2f' ...            % 3: ppt
%             '\t%g\t%g\t%g\t%g' ... % 4-7
%             '\t%g\t%d\t%g\t%d' ... % 8-11
%             '\t%d\t%d\t%d\t%g' ... % 12-15
%             '\t%g\t%g\t%d\t%g' ... % 16-19
%             '\t%d\t%g\t%g\t%g' ... % 20-23
%             '\t%g\t%g\t%g\t%g' ... % 24-27
%             '\t%d\t%g\t%g\t%g\n']; % 28-31
%     otherwise
%         error('Sitename not matched');
% end

OUTFILE=[sitename{1} '_MET_' num2str(jd_start) '_' num2str(jd_end-1) '_' num2str(met.startyear) '.dat'];

%% [2] select row range
%range_of_interest=intersect( find(met.decday <= jd_end), ...
%    find(met.decday > jd_start) );

% generate timebase
[ce_year ce_day ce_hrmin ce_hr1 ce_min1]=c_eur_gen_timebase(met);
% check for missing values / time base errors.
[timebase_status]=c_eur_check_timebase(ce_day, ce_hr1, ce_min1);
if (timebase_status > 0)
    error (['Timebase error. Stopping.']);
end


%% [4] column swapping
% desired order for CARBO EUROPE is:
%
% -1- precipitation
%global or short wave incoming radiation	%
%reflected or short wave outgoing radiation
%long wave incoming radiation
%long wave outgoing radiation
%net radiation
%diffuse radiation (global)
% -8- photosynthetic photon flux density
%diffuse photosynthetic photon flux density
%reflected photosynthetic photon flux density
%below-canopy photosynthetic photon flux density
%light interception
% -13- air temperature
%pressure
%canopy radiative temperature
%bole temperature
%soil temperature superficial
%soil temperature medium
%soil temperature deep
%soil water content superficial
%soil water content medium
%soil water content deep
%soil heat flux g1
% soil heat flux	g2
% -25- relative humidity
%Wind direction
%Wind horizontal speed
%ADDITIONAL PARAMETER	[soil temp 2]
%ADDITIONAL PARAMETER	[swc 15cm]
%ADDITIONAL PARAMETER
%ADDITIONAL PARAMETER
%ADDITIONAL PARAMETER
%ADDITIONAL PARAMETER	%
%ADDITIONAL PARAMETER


%% [3] manipulation.

% pressure in kpa :
met.Pkpa.data=met.P.data./10;
met.Pkpa.unit='kPa';

empty = ones(1,numel(met.ppt.data)).*MISSING_VALUE;   % array of -9999s for values not measured.

% all SWC in % vol (not fractions) :
SWC_all_fields={'Theta_5','Theta_15','Theta_25','Theta_20','Theta_30',...
    'Theta_40','Theta_50','Theta_0_30','Tsoil_75mm'};
for i_field=1:numel(SWC_all_fields)
    if (isfield(met,SWC_all_fields{i_field}))
        curr_field=eval(['met.',SWC_all_fields{i_field}]);
        if (isfield(curr_field,'data'));
            evalstr=['met.',SWC_all_fields{i_field},'.data=met.',SWC_all_fields{i_field},'.data.*100'];
        
            eval(evalstr);
        end
        evalstr=['met.',SWC_all_fields{i_field},'.unit=''%'''];
        eval(evalstr);
     
    else
        % insert a dummy field if not found:
        eval(['met.',SWC_all_fields{i_field},'.data=empty.''']);
    end
end

% RH in %
met.RH.data=100.*met.RH.data;

% correct G for soil heat storage (is this necessary for C/EUR?)
%met.G1.data=storage_soil_heat_flux(met.G1.data,met.Tsoil1.data,met.Theta_5.data,0.05,0.5);
%met.G2.data=storage_soil_heat_flux(met.G2.data,met.Tsoil2.data,met.Theta_5.data,0.05,0.5);

% make the Rswin, Rlwin, PAR baseline values zero
[par_corrected zero_error pc_rejected]=correct_rad_zero_error...
    (met.PAR.data,met.decday,1);
met.PAR.data=par_corrected;
[rswin_corrected zero_error pc_rejected]=correct_rad_zero_error...
    (met.Rswin.data,met.decday,5);
met.Rswin.data=rswin_corrected;

% apply some of the carbo europe filters (radiation, extrme values etc)
limits=define_met_carbo_europe_limits();
met=met_carbo_europe_filters(met,limits);


met_out = [met.ppt.data.';     met.Rswin.data.' ; met.Rswout.data.'  ;  met.Rlwin.data.'; met.Rlwout.data.';    met.Rn.data.';    empty; met.PAR.data.'; ...
    empty; empty; empty; empty;     met.T_gauge.data.' ;   met.Pkpa.data.'   ;  met.T_surf.data.';   empty; ...
    met.Tsoil1.data.' ;   empty; met.Tsoil_75mm.data.' ;    met.Theta_5.data.' ;  met.Theta_25.data.' ; met.Theta_50.data.' ;...
    met.G1.data.' ; met.G2.data.' ;   met.RH.data.';    met.propwinddir.data.'; met.propwindspeed.data.';   met.Tsoil2.data.'; met.Theta_15.data.'].';
% %   23      24      25    26        27          28      29

% prepend the julian day & hours info:
met_out=horzcat(reshape(ce_hr1,numel(ce_hr1),1), reshape(ce_min1,numel(ce_min1),1), met_out);

% now replace all NaNs with -9999 for carbo europe
met_out(find(isnan(met_out)))=MISSING_VALUE;
% visual check: call plot_met_out.m to plot all the met params.
% plot_met_out(met_out);

%% [6] write to output file. append to header file.
format_ce = '%10s\t%02d:%02d\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n';

status=copyfile(header_file,OUTFILE);
fid=fopen(OUTFILE,'a');
for i = 1:size(met_out,1)
    datestr = jday2ddmmyyyy(ce_year(i),ce_day(i));
    % format_str is site dpenedent. defined in switch statement @ start
    fprintf(fid,format_ce,datestr,met_out(i,:));
end
fclose(fid);
disp(['Finished. ',num2str(size(met_out,1)),' values written to ' OUTFILE]);
result =0;