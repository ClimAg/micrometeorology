% define_met_carbo_europe_limits.m
% PL 26.08.2006
%
% Set up a limits structure for CARBOEUROPE met data filtering.
% See http://gaia.agraria.unitus.it/ or template s/sheet for allowed limits.
%
% NB Many limits are still undefined!

function limits = define_met_carbo_europe_limits()

limits.T_surf  = [-40 50];      % deg C
limits.Rn = [-200 1000];    % W m-2
limits.G = [-300 300]; % W m-2

limits.Rswin= [-50 1200]; % W m-2
limits.Rswout=[-50 1200];

limits.T_gauge=[-40 50]; % deg C

limits.Pkpa=[70 130]; % kPa
limits.T_soil=[-30 40]; % deg C
limits.RH=[5 100]; % percent
limits.SWC=[0 70]; % percent vol
