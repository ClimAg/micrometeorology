% met_carbo_europe_filters.m
% PL 06.01.2006
% Apply some filters for CARBO EUROPE met data.
% The filters are defined in the Excel spreadsheet template met_cip.xls
%    and in the notes on grassland stations.
% PL 26.08.2006 Updated to filter Rswin/out, RH, P, Tsoil etc.
% PL 21.11.2007 Added check for existence of Pkpa field in met struc.t
%
% INPUTS:
%   met:    a met data structure from, e.g., read_cork_data
%   limits: a structure containing the limits for met quantities
%               (see define_met_carbo_europe_limits.m )
%   MISSING_VALUE: the value to replace bad values with (usually -9999)
%
function met = met_carbo_europe_filters(met, limits, MISSING_VALUE)

if (~exist('MISSING_VALUE', 'var'))
    MISSING_VALUE = -9999;
end

% Rn filter : -200 < Rn < 1000 W /m2
met.Rn.data(find(met.Rn.data<limits.Rn(1))) = MISSING_VALUE;
met.Rn.data(find(met.Rn.data>limits.Rn(2))) = MISSING_VALUE;

% Rsw in/out filter -50, 1200. Known as Rg, Rr in Carboeurope.
met.Rswin.data(find(met.Rswin.data<limits.Rswin(1))) = MISSING_VALUE;
met.Rswin.data(find(met.Rswin.data>limits.Rswin(2))) = MISSING_VALUE;
met.Rswout.data(find(met.Rswout.data<limits.Rswout(1))) = MISSING_VALUE;
met.Rswout.data(find(met.Rswout.data>limits.Rswout(2))) = MISSING_VALUE;


% mark values where Rswout > Rswin as bad (allow for small zero offset in
% both)
bad_sw = find(met.Rswout.data > met.Rswin.data);
non_zero = find(met.Rswout.data > 10);
exclude_sw = intersect(non_zero,bad_sw);
met.Rswout.data(exclude_sw) = NaN;
met.Rswin.data(exclude_sw)=NaN;

% mark values where Rlwin > Rlwout as bad
bad_lw = find (met.Rlwin.data > met.Rlwout.data);
met.Rlwout.data(bad_lw) = NaN;
met.Rlwin.data(bad_lw) = NaN;

% T_surf / canopy temperature:
met.T_surf.data(find(met.T_surf.data<limits.T_surf(1))) = MISSING_VALUE;
met.T_surf.data(find(met.T_surf.data>limits.T_surf(2))) = MISSING_VALUE;

% air temperature
met.T_gauge.data(find(met.T_gauge.data<limits.T_gauge(1))) = MISSING_VALUE;
met.T_gauge.data(find(met.T_gauge.data>limits.T_gauge(2))) = MISSING_VALUE;

% G (soil heat flux)
met.G1.data(find(met.G1.data<limits.G(1))) = MISSING_VALUE;
met.G2.data(find(met.G2.data<limits.G(1))) = MISSING_VALUE;
met.G1.data(find(met.G1.data>limits.G(2))) = MISSING_VALUE;
met.G2.data(find(met.G2.data>limits.G(2))) = MISSING_VALUE;

% Pressure limits (kPa) if Pkpa field is present
if (isfield(met,'Pkpa'))
    met.Pkpa.data(find(met.Pkpa.data<limits.Pkpa(1))) = MISSING_VALUE;
    met.Pkpa.data(find(met.Pkpa.data>limits.Pkpa(2))) = MISSING_VALUE;
end
% Soil temperatures
% Tsoil 75 mm (if the field exists)
if (isfield(met,'Tsoil_75mm'))
    met.Tsoil_75mm.data(find(met.Tsoil_75mm.data<limits.T_soil(1))) = MISSING_VALUE;
    met.Tsoil_75mm.data(find(met.Tsoil_75mm.data>limits.T_soil(2))) = MISSING_VALUE;
end

met.Tsoil1.data(find(met.Tsoil1.data<limits.T_soil(1))) = MISSING_VALUE;
met.Tsoil1.data(find(met.Tsoil1.data>limits.T_soil(2))) = MISSING_VALUE;

met.Tsoil2.data(find(met.Tsoil2.data<limits.T_soil(1))) = MISSING_VALUE;
met.Tsoil2.data(find(met.Tsoil2.data>limits.T_soil(2))) = MISSING_VALUE;

% RH
met.RH.data(find(met.RH.data<limits.RH(1))) = MISSING_VALUE;
met.RH.data(find(met.RH.data>limits.RH(2))) = MISSING_VALUE;

% soil moisture content
met.Theta_5.data(find(met.Theta_5.data<limits.SWC(1))) = MISSING_VALUE;
met.Theta_5.data(find(met.Theta_5.data>limits.SWC(2))) = MISSING_VALUE;

if (isfield(met,'Theta_15'))
    met.Theta_15.data(find(met.Theta_15.data<limits.SWC(1))) = MISSING_VALUE;
    met.Theta_15.data(find(met.Theta_15.data>limits.SWC(2))) = MISSING_VALUE;
end

if (isfield(met,'Theta_25'))
    met.Theta_25.data(find(met.Theta_25.data<limits.SWC(1))) = MISSING_VALUE;
    met.Theta_25.data(find(met.Theta_25.data>limits.SWC(2))) = MISSING_VALUE;
end

met.Theta_50.data(find(met.Theta_50.data<limits.SWC(1))) = MISSING_VALUE;
met.Theta_50.data(find(met.Theta_50.data>limits.SWC(2))) = MISSING_VALUE;

if (isfield(limits,'ppt'))
    met.ppt.data(find(met.ppt.data<limits.ppt(1))) = MISSING_VALUE;
    met.ppt.data(find(met.ppt.data>limits.ppt(2))) = MISSING_VALUE;
end
if ((isfield(limits,'ppt')) && (isfield(met,'ppt2')))
    met.ppt2.data(find(met.ppt2.data<limits.ppt(1))) = MISSING_VALUE;
    met.ppt2.data(find(met.ppt2.data>limits.ppt(2))) = MISSING_VALUE;
end
if ((isfield(limits,'PAR')) && (isfield(met,'PAR')))
    met.PAR.data(find(met.PAR.data<limits.PAR(1))) = MISSING_VALUE;
    met.PAR.data(find(met.ppt.data>limits.PAR(2))) = MISSING_VALUE;
end
