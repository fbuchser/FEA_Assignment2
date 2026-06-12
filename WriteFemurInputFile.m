% Function call:
% WriteFemurInputFile(['data/material_mapped_mesh.txt'], 'data/','B')

function WriteFemurInputFile(GEMfilename, OutputPath, Method)

if exist(GEMfilename,'file')==0
    error('*****Can not find input file %s*****', GEMfilename)
end

if exist(OutputPath,'dir')==0
    error('*****Can not find output path %s*****', OutputPath)
end

%----------- BEGIN READING INPUT FILE FROM MITK-GEM ----------
fprintf('*\n')
fprintf('*****READING ANSYS INPUT FILE*****\n')
fprintf('*\n')
[~, FileName, ext] = fileparts(GEMfilename);
FileName = [FileName, ext];

% ******************OPEN FILE*****************
fprintf('*****INPUT FILE: %s*****\n', FileName)
fprintf('*\n')
fid0=fopen(GEMfilename);

% ******************READ NODES*****************
tline='!!!';

fprintf('*****READING NODE DATA*****\n')
fprintf('*\n')
while length(tline)<6 || ~strcmpi(tline(1:12),'#BEGIN NODES')
    tline=fgetl(fid0);
end

fgetl(fid0); % Reading 1 more header line
fgetl(fid0); % Reading 1 more header line

NODES = cell2mat(textscan(fid0, '%f,%f,%f,%f,%f', 'EndOfLine','\r\n'));
fgetl(fid0); % Reading '#END NODES'

% ******************READ ELEMENTS*****************
tline='!!!';

fprintf('*****READING ELEMENT DATA*****\n')
fprintf('*\n')
while length(tline)<6 || ~strcmpi(tline(1:15),'#BEGIN ELEMENTS')
    tline=fgetl(fid0);
end

fgetl(fid0); % Reading 1 more header line
fgetl(fid0); % Reading 1 more header line

ELEMENTS = cell2mat(textscan(fid0, '%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f', 'EndOfLine','\r\n'));
fgetl(fid0); % Reading '#END ELEMENTS'

% ******************BUILD nodes AND elements*****************
nyy = 0.3;
nodes = NODES(:,1:4);

% Method A: node-averaged E
% Method B: element-averaged E
switch Method
    case 'A'
        %Write the Method A material file
        El_E = ELEMENTS(:,12);
    case 'B'
        %Write the Method B material file
        El_E = ELEMENTS(:,13);
end

elements = [ELEMENTS(:,1:5) El_E nyy*ones(length(El_E(:,1)),1)];

% Keep only nodes referenced by elements
el_col = elements(:,2:5);
el_col = el_col(:);
I = intersect(nodes(:,1),el_col);
nodes = nodes(I,:);

save(fullfile(OutputPath, 'FemurInput.mat'), 'nodes','elements')
    