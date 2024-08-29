%c_eur_gen_timebase.m
% PL 21.12.2005
%
% Generate the year , day, hrmin, hr, min vectors from a flux structure.
% For use with carbo europe flux/met files.
%
% INPUTS: (f) flux or met structure from postproc;
%           [optional:] i_range: list of indices of values to consider.
% OUTPUTS: [ce_year ce_day ce_hrmin ce_hr1 ce_min1]


% usage: [ce_year ce_day ce_hrmin ce_hr1 ce_min1]=c_eur_gen_timebase(f,i_range);

function [ce_year ce_day ce_hrmin ce_hr1 ce_min1] = c_eur_gen_timebase(f, i_range)
if (~exist('i_range'))
    i_range=1:numel(f.decday);
end
ce_year = zeros(numel(i_range),1)+f.startyear;
ce_day = floor(f.decday(i_range));
ce_dechr = (0.5.*round(48.*(f.decday(i_range)-ce_day)));
ce_hr = floor(ce_dechr);
ce_min = 60.*(ce_dechr-ce_hr);
ce_hrmin = ce_hr.*100+ce_min; % hhmm format

%change 0000 to 2400 if necessary
ce_hrmin(find(ce_hrmin==0))=2400;
ce_hr1=floor(ce_hrmin/100);       % round time1 to the nearest integer < or equal to time1
ce_day(find(ce_hrmin==2400))=ce_day(find(ce_hrmin==2400)) -1; % time 2400 "belongs" to prev. day.
%ce_min1=ce_hrmin-hr*100;         % determine the minutes

ce_min1=ce_hrmin-ce_hr1.*100;         % determine the minutes
