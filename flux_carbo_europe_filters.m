% flux_carbo_europe_filters.m
% PL 27.07.2004
% modified 09.08.2004 to perform size check on input data.
% modified 22.12.2005   : (1) "stab" -> "stability" 
%                       : (2) filters NaNs (moved here from flux_carbo_europe_2005.m) 
%                       : (3) updated structures to reflect new
%                       nomenclature
% PL 12.11.2007 : Changed calc_itc() call syntax to new version
% PL 13.02.2009 : Initial f.wc.qc now set to 1 (formerly 2)
%
% Applies the  data filters detailed in the document "instruction on data delivery and database" from Sabina Dore.
% 
% these filters replace values outside the allowable ranges with -9999. 
%
% INPUTS: a structure with the following fields
% f.Cbar.data               CO2 conc
% f.wc.data          CO2 flux
% f.wT.data            sensible heat
% f.qT.data           latent heat flux
% f.ustar.data       frictional velocity m/s
% f.qbar.data           water vapour conc mmol/mol
% f.stability.data        stability paramtetr  z/L
%
% OUTPUTS: the same structure with values replaced by -9999 if outside the
% ranges.
% new field added : "f.wc.qc" indicates quailty control flag
% also "f.H.qc", "f.LE.qc"

function f = flux_carbo_europe_filters(f, limits, MISSING_VALUE)
% initialisateions:
if (~exist('limits'))
    % define the filter limits:
    limits.Cbar  = [300 700];      % ppm (umol/mol)
    limits.wc   = [-50  50];
    limits.H    = [-200 800];
    limits.LE   = [-200 800];
    limits.ustar= [-1   3.5];   % m/s
    limits.qbar    = [0.001 40];   % mmol/mol (ppm_mol)
    limits.stability = [-50 50];
    limits.tau       = [0 5]; % kg m-1 s-2
    limits.ustarthreshold = 0.2; %m/s
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (~exist('MISSING_VALUE'))
    MISSING_VALUE = -9999;
end

if (...
        (numel(f.Cbar.data) ~= numel(f.wc.data)) || (numel(f.Cbar.data) ~= numel(f.H.data)) ...
        || (numel(f.Cbar.data) ~= numel(f.LE.data)) || (numel(f.Cbar.data) ~= numel(f.ustar.data)) || (numel(f.Cbar.data) ~= numel(f.qbar.data)) ...
        || (numel(f.Cbar.data) ~= numel(f.stability.data))...
        )
    error('Dimensions of input variables are not equal.');
end
% stability filtering:
f.stability.data(find(f.stability.data < limits.stability(1))) = NaN;
f.stability.data(find(f.stability.data > limits.stability(2))) = NaN;


% carry out the ITC test: (Foken instructions)
[ITC meas_itc mod_itc]= calc_itc(f.stability.data,...
    f.ustd.data,...
    f.vstd.data,...
    f.wstd.data,...
    f.Tstd.data,...
    f.ustar.data,...
    f.wT.data);
%ITC = calc_itc(f.stability.data, f.wstd.data, f.ustar.data);
f.stability.data(find(isnan(f.stability.data))) = MISSING_VALUE;

f.wc.qc = ones(numel(ITC),1)*1;    % initialise all values to one ("questionalbe" quality)
f.wc.qc(find(ITC >= 100)) = 2;
f.wc.qc(find(ITC < 100)) = 1;      % order is important herE!
f.wc.qc(find(ITC < 30)) = 0;

f.tau.qc=f.wc.qc;  % tau shares the same 'bad' values as wc so far.

% u* filtering: (Vesna thesis) (threhsold u*  = 0.2 ms-1)
umask=(f.ustar.data<limits.ustarthreshold);     
umask=umask.*2;  % must be >=2 to activate Reichstein gap filler
f.wc.qc = f.wc.qc+umask; % combine with existing qc value 
% max allowed quality flag value for carbo europe is 2 :
f.wc.qc(find(f.wc.qc > 2)) = 2;
% all the above tests also apply to H, LE:
f.H.qc = f.wc.qc;
f.LE.qc = f.wc.qc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% perform the abs value filtering:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Cbar, CO2 conc: 
f.Cbar.data(find(f.Cbar.data < limits.Cbar(1))) = MISSING_VALUE;
f.Cbar.data(find(f.Cbar.data > limits.Cbar(2))) = MISSING_VALUE;
f.Cbar.data(find(isnan(f.Cbar.data))) = MISSING_VALUE;

% wc, CO2 flux:
f.wc.data(find(f.wc.data < limits.wc(1))) = MISSING_VALUE;
f.wc.data(find(f.wc.data > limits.wc(2))) = MISSING_VALUE;
f.wc.data(find(isnan(f.wc.data))) = MISSING_VALUE;
f.wc.qc(find(f.wc.data==MISSING_VALUE)) = MISSING_VALUE; % if wc is bad, set qc flag to 2

% fc_gf, gapfilled CO2 flux:
% f.fc_gf(find(f.fc_gf < limits.fc_gf(1))) = MISSING_VALUE;
% f.fc_gf(find(f.fc_gf > limits.fc_gf(2))) = MISSING_VALUE;
% f.fc_gf(find(isnan(f.fc_gf))) = MISSING_VALUE;

% tau (uw), momentum flux:
f.tau.data(find(f.tau.data < limits.tau(1))) =  MISSING_VALUE;
f.tau.data(find(f.tau.data > limits.tau(2))) =  MISSING_VALUE;
f.tau.data(find(isnan(f.tau.data))) = MISSING_VALUE;
f.tau.qc(find(f.tau.data == MISSING_VALUE)) = MISSING_VALUE;

% H, sensible heat flux
f.H.data (find(f.H.data < limits.H(1))) = MISSING_VALUE;
f.H.data (find(f.H.data  > limits.H(2))) = MISSING_VALUE;
f.H.data (find(isnan(f.H.data ))) = MISSING_VALUE;
f.H.qc(find(f.H.data==MISSING_VALUE)) = MISSING_VALUE; 
%( if fluxes are -9999 the quality flags should also be -9999 (Sabina Dore))

% LE, latent heat flux:
f.LE.data (find(f.LE.data  < limits.LE(1))) = MISSING_VALUE;
f.LE.data (find(f.LE.data  > limits.LE(2))) = MISSING_VALUE;
f.LE.data (find(isnan(f.LE.data))) = MISSING_VALUE;
f.LE.qc(find(f.LE.data==MISSING_VALUE)) = MISSING_VALUE; 

% u*, fricitional velocyti:
f.ustar.data (find(f.ustar.data  < limits.ustar(1))) = MISSING_VALUE;
f.ustar.data (find(f.ustar.data  > limits.ustar(2))) = MISSING_VALUE;
f.ustar.data (find(isnan(f.ustar.data))) = MISSING_VALUE;

% qbar, H2O concentration:
f.qbar.data(find(f.qbar.data < limits.qbar(1))) = MISSING_VALUE;
f.qbar.data(find(f.qbar.data > limits.qbar(2))) = MISSING_VALUE;
f.qbar.data(find(isnan(f.qbar.data))) = MISSING_VALUE;
