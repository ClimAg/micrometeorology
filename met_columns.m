% met_columns.m
% PL 02.07.2004
% Produces the met data column mapping for Dripsey/Wexford/Kerry.
%
% inputs:   sitename '1=D','2=K','3=W'
% returns: a variable with column indices to map met data to carbo europe
% format data columns
function column_map = met_columns(sitename)

switch sitename
    case 1  % Dripsey
        column_map = [5 7 9 8 10 6 -1 28 -1 -1 ...
                -1 -1 12 21 15 -1 18 -1 20 24 ...
                26 27 16 17 11 -1 31 19 25];
   
    case 2  % Kerry
        column_map = [5 8 10 9 11 7 -1 12 -1 -1 ...
                -1 -1 14 24 15 -1 18 20 23 25 ...
                27 30 16 17 13 34 33]; % same as W except last 2 vals
    case 3  % WExford
        
        column_map = [5 8 10 9 11 7 -1 12 -1 -1 ...
                -1 -1 14 24 15 -1 18 20 23 25 ...
                27 30 16 17 13 -1 -1];
    otherwise
        column_map = [];
end
% desired order is:
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
%bole temperature (tree trunk)	
%soil temperature superficial 	
%soil temperature medium	
%soil temperature deep	
% -20- soil water content superficial	
%soil water content medium	
%soil water content deep	
%soil heat flux g1	
% soil heat flux	g2
% -25- relative humidity 	
%Wind direction 	
%Wind horizontal speed 	
%ADDITIONAL PARAMETER	
%ADDITIONAL PARAMETER	
%ADDITIONAL PARAMETER	
%ADDITIONAL PARAMETER	
%ADDITIONAL PARAMETER	
%ADDITIONAL PARAMETER	%
%ADDITIONAL PARAMETER
