%flux_carbo_europe.m
% process the flux data from Dripsey/Wexford/Kerry for Carbo Europe
% PL 1/7/2004
% PL 11-12.2005 This version modified to use the new flux data structures.
% from Slava's new postproc file.
% PL 26.05.2006 Use 0 for storage terms, not 9999 (D. Papale email, May'06)
% PL 12.11.2007 Changed call to calc_itc (new syntax). Added 2007 code.
% PL 22.11.2007 Use met T, P when available.
%               Changed default molar vol from 22.4 to 23.9 L.
% PL 27.11.2007 Added code to handle low goodpts intervals in 2005/06.
% PL 13.09.2009 Renamed to flux_carbo_europe.m; converted to function.
%
% [1]   read in the postprocessed data.
% [2]   select a time range
% [3]   manipluation:
%    change 0000 to 2400 if necessary
% [4]   swap columns (to get correct carbo europe order
% [5]   do checking
% [6]   write to output file ready for carbo europe
%
% INPUTS
% f     flux structure, postprocessed, filtered. 
% F     flux structure, derived values
% met   met structure
% s     site structure
% ignore_low_goodpts  set to 1 to remove goodpts filtering
%
function result=flux_carbo_europe(f,F,met,s,ignore_low_goodpts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  general settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load_constants; % slava's physical constants (for gas const.)
startup; % path names
sitename = s.name; % 
% define time range of interest
        jd_start = 1;   % from AND? including
        %jd_end = 182;   % up to AND including
        jd_end = 366;
%ROOT_DIR = 'E:\Data\CARBOEUROPE\';
header_file=fullfile(PATH_CARBOEUROPE,'flux_header.txt');

MISSING_VALUE=-9999;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  end of general settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% switch sitename
%     case 'D'
%       
%         num_cols=29;
%         full_sitename='IE_DRI';
%     otherwise
%         error('Sitename not matched');
% end
percent_nan=100.*sum(isnan(f.wc.data)./numel(f.wc.data));
disp([num2str(percent_nan),' % NaN after postproc.']);

% [2] select row range
range_of_interest=intersect( find(f.decday <= jd_end), ...
    find(f.decday > jd_start) );


% generate timebase 
[ce_year ce_day ce_hrmin ce_hr1 ce_min1]=c_eur_gen_timebase(f);
% check for missing values / time base errors.
[timebase_status]=c_eur_check_timebase(ce_day, ce_hr1, ce_min1);
if (timebase_status > 0)
    error (['Timebase error. Stopping.']);
end

% use gapfilled data if available; otherwise use NaN
if (~isfield(f.wc,'gapfilled'))
    f.wc.gapfilled.data=ones(numel(f.decday),1).*MISSING_VALUE;
    f.wc.gapfilled.unit='micromol m-2 s-1';

end

% add the stability parameter and H, LE to the "f" structure (for use in filters):
f.stability.data=F.stability.data;
f.H=F.heat;
f.LE=F.evaporation;
f.ustar=F.ustar;
% calculate the vpd:
P_default=991; % mbar
molar_vol_default=23.9e-3; % [L] (based on 991 mbar, 285 K average Dripsey conditions)
% replace NaNs with -9999 amd change units from mmol m-2 s-1 to umol m-2 s-1 before filtering:
f.wc.data=errs_to_nan(f.wc.data);   % replace all error codes (9999 etc) with NaN
f.wc.data=1000.*f.wc.data;
f.wc.data(isnan(f.wc.data))= MISSING_VALUE; % now replace all NaNs with -9999
% this seems strange but first operation is many-to-one mapping, second op
% is one-to-one. this is a roundabout way of converting all errs to -9999.

press=P_default; % default atm pressure, defined at top of funcution
if (isfield(met,'P'))
    % replace default value with measured P
    press=ones(size(f.wc.data)).*P_default;
    % find overlapping flux/met data
    met_within_range=intersect(find(met.datenum>=f.datenum(1)), ...
             find(met.datenum<=f.datenum(end)));
    % filter absolute P values outside [900,1100]
    met.P.data(find((abs(met.P.data-1000))>100))=NaN;
    % filter absolute Tgauge values outside [-5 35]
    met.T_gauge.data(find((abs(met.T_gauge.data-15))>20))=NaN;
    % resample the data
    T_i = interp1(met.datenum,met.T_gauge.data,f.datenum,'linear');
    press_i = interp1(met.datenum,met.P.data,f.datenum,'linear'); 
    press=1e2.*press_i; % convert from millibar to Pa
    T=T_i+273.15; % convert from deg C to K
    molar_vol=R.*T./press;
    % PL comment : should use a moving average, not overall avg! 
    bad_molar_vol=(find(isnan(molar_vol)));
    molar_vol(bad_molar_vol)=nanmean(molar_vol);
end
molar_vol=reshape(molar_vol,size(f.Cbar.data)); % units : m^3
% change units from umol m-3 to ppm_v before filtering:
% [PL 26.11.2007 : units of Cbar appear to be mmol m-3 beforehand]
f.Cbar.data=molar_vol.*f.Cbar.data.*1000;
% change units of water vapor here too (from mmol m-3 to mmol mol-1)
f.qbar.data=molar_vol.*f.qbar.data;

% add the momentum flux to the "f" structure % tau (momentum flux in kg m-1 s-2)
rho_a=1.22; % kg m-3
f.tau.data=-1.*rho_a.*f.uw.data; % PL : seems to be a sign convention mismatch betw Slava and C/Eur
f.tau.unit='kg m-1 s-2';

%  filter out values outside acceptable range for fc, C , H LE , ustar,
%  fc_gf:
%  (replacing with -9999)
% warning: this will possibly re-introduce gaps into the filled data! 
% To prevent this from happening, go back and filter the data BEFORE it is gap filled. 
f = flux_carbo_europe_filters(f);

percent_9999=100.*sum((f.wc.data==MISSING_VALUE)./numel(f.wc.data));
disp([num2str(percent_9999),' % -9999 after c/eur filters.']);

% PL extra filter required during datalogger malfunction 06/07:
% set qc >= 1 when goodpts was below 17000. 
if (ignore_low_goodpts==1)
    failed_goodpts=(f.goodpts<17000);
    f.wc.qc     =(failed_goodpts.*f.wc.qc); % multiply to preserve existing vlaues
    f.LE.qc     =(failed_goodpts.*f.LE.qc);
    f.H.qc      =(failed_goodpts.*f.H.qc);
    f.tau.qc    =(failed_goodpts.*f.tau.qc);
end

% UCC filters (ppt, extreme values): [Rswin threshold =25, no ppt filtering]
% [good_index bad_index]=co2_flux_filter(f,met,site,25,9999);
% good_wc=find(f.wc.qc == 0);
% f.wc.qc(intersect(bad_index,good_wc))=1;

%% [6] write to output file.
empty = ones(numel(f.Cbar.data),1).*MISSING_VALUE;   % array of -9999s for values not measured.
zero_list=zeros(numel(f.Cbar.data),1);

OUTFILE=[sitename{1},'_FLUX_', num2str(jd_start),'_',num2str(jd_end-1),'_',num2str(f.startyear),'.dat'];
% [6] write to output file. append to header file.
status=copyfile(header_file,OUTFILE);
fid=fopen(OUTFILE,'w');
format_ce = '%10s\t%02d:%02d\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n';
for i_range = 1:numel(range_of_interest)
    i=range_of_interest(i_range); % current timestep 
    datestr = jday2ddmmyyyy(ce_year(i),ce_day(i));
    fprintf(fid,format_ce,datestr,ce_hr1(i),ce_min1(i),...
        f.Cbar.data(i),...                 % col 3
        f.qbar.data(i),...                 % col 4
        f.stability.data(i),...              % etc...
        f.wc.data(i),...
        f.wc.qc(i),...
        f.H.data(i),...
        f.H.qc(i),...
        f.LE.data(i),...                    % 10
        f.LE.qc(i),...
        f.tau.data(i),...
        f.tau.qc(i),...
        f.ustar.data(i),...
        empty(i),...
        zero_list(i),... % canopy storage term
        zero_list(i),... % 17 canopy storage term
        zero_list(i),... % canopy storage term
        zero_list(i),... % canopy storage term
        f.gapfilled.wc.data(i)...            % 20 gapfilled co2 flux
        );
end
fclose(fid);


disp(['Finished. ',num2str(numel(range_of_interest)),' values written to ' OUTFILE]);

% APPENDIX : CARBO EUROPE desired data order:
%
%-1-date,
%time,
%carbon dioxide concentration ,
%Water vapour concentration ,
%-5-atmosphere stability parameter ,
%carbon dioxide ,
%quality class,
%sensible heat ,
%quality class,
%-10-latent heat ,
%-11-quality class,
%Momentum,
%quality class,
%friction velocity,
%tree transpiration ,
%-16-canopy heat storage ,
%CO2 storage in canopy air layer ,
%LATENT HEAT in canopy air layer ,
%heat storage in canopy air layer ,
%-20- gapfilled flux data ,
%-21-ADDITIONAL PARAMETER,
%ADDITIONAL PARAMETER,
%ADDITIONAL PARAMETER,
%ADDITIONAL PARAMETER,
%ADDITIONAL PARAMETER

% postproc: Outarray_carbo_europe = [year day hrmin C q stab fc qc_gf H qc.' LE qc_gf empty empty u_star empty empty empty empty empty empty vpd fc_orig];
result=0;