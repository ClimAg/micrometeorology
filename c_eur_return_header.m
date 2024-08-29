% c_eur_return_header.m
% PL 22.12.2005
% Return a carbo-europe format file header, either MET or FLUX
%%
% INPUTS: type is either 'MET' or 'FLUX'
%
function headerstring = c_eur_return_header(type)
type=upper(type);
switch type
    case 'FLUX'
        headerfilename='E:\Data\CARBOEUROPE\documents_templates\flux_header.txt';
    case 'MET'
        headerfilename='E:\Data\CARBOEUROPE\documents_templates\met_header.txt';
    otherwise
        error('Unrecognised argument. Use either FLUX or MET');
end
fid=fopen(headerfilename,'r');
headerstring1 = fgetl(fid);
headerstring2 = fgetl(fid);
headerstring3 = fgetl(fid);
headerstring = strvcat(headerstring1, headerstring2,headerstring3);
fclose(fid);
