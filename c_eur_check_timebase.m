% c_eur_check_timebase.m
% PL check the flux/met timebase for errors.
% Must adhere to the carbo europe standard.
%
% INPUTS: (ce_day, ce_hr1, ce_min1) - all from c_eur_gen_timebase()
%
% OUTPUTS: timebase_status - 1=>bad ; 0=>good.
%           bad_points : list 1/0 indicating bad points.
%
% Usage: [timebase_status]=c_eur_check_timebase(ce_day, ce_hr1, ce_min1);
function [timebase_status bad_points] = c_eur_check_timebase(ce_day, ce_hr1, ce_min1)
timebase_status=0; % initialise.
logged_timebase = ce_day+ (ce_hr1 + ce_min1./60)./24; 
calc_timebase = logged_timebase(1) + linspace(0,(numel(ce_day)-1)./48,...
        numel(ce_day));
    calc_timebase=reshape(calc_timebase,size(logged_timebase));
bad_points = (abs(logged_timebase-calc_timebase)> 0.002);   % 0.002 is 0.1 * time increment.
if (sum(bad_points) > 0)
    timebase_status=1;
end
