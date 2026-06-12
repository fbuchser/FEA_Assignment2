%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%WriteFemurInputFile_rev2(['material mapped mesh 6mm.txt'],'B')
%WriteFemurInputFile_rev2(['material mapped mesh 3mm.txt'],'B')
function WriteFemurInputFile_rev2(GEMfilename,Method)

if exist(GEMfilename,'file')==0
    error('*****Can not find input file %s*****',GEMfilename)
end


%----------- BEGIN READING INPUT FILE FROM MITK-GEM ----------
%----------- BEGIN READING INPUT FILE FROM MITK-GEM ----------
fprintf('*\n')
fprintf('*****READING ANSYS INPUT FILE*****\n')
fprintf('*\n')
[PathName,FileName,ext]=fileparts(GEMfilename);
FileName = [FileName,ext];

if isnumeric(FileName)
    error('*****Operation Cancelled*****\n\n')
end

% ******************OPEN FILE*****************
fprintf('*****INPUT FILE: %s*****\n',FileName)
fprintf('*\n')
fid0=fopen(GEMfilename);


% ******************READ NODES*****************
tline='!!!';
NODES = [];
fprintf('*****READING NODE DATA*****\n')
fprintf('*\n')
while length(tline)<6 || ~strcmpi(tline(1:12),'#BEGIN NODES')
    tline=fgetl(fid0);
end
tline=fgetl(fid0); % Reading 1 more header line
tline=fgetl(fid0); % Reading 1 more header line
NODES = cell2mat(textscan(fid0, '%f,%f,%f,%f,%f', 'EndOfLine','\r\n'));
fgetl(fid0); % Reading '#END NODES'


% ******************READ ELEMENTS*****************
tline='!!!';
ELEMENTS = [];
fprintf('*****READING ELEMENT DATA*****\n')
fprintf('*\n')
while length(tline)<6 || ~strcmpi(tline(1:15),'#BEGIN ELEMENTS')
    tline=fgetl(fid0);
end
tline=fgetl(fid0); % Reading 1 more header line
tline=fgetl(fid0); % Reading 1 more header line
ELEMENTS = cell2mat(textscan(fid0, '%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f', 'EndOfLine','\r\n'));
fgetl(fid0); % Reading '#END ELEMENTS'
%index = [1 2 3 4 5 6 7 8 11 9 10 12 13];
%ELEMENTS = ELEMENTS(:,index);


% ******************READ SURFACES*****************
tline='!!!';
SURFACES = [];
fprintf('*****READING EXTERNAL SURFACE DATA*****\n')
fprintf('*\n')
while length(tline)<6 || ~strcmpi(tline(1:14),'#BEGIN SURFACE')
    tline=fgetl(fid0);
end
tline=fgetl(fid0); % Reading 1 more header line
SURFACES = cell2mat(textscan(fid0, '%f,%f,%f,%f,%f,%f,%f', 'EndOfLine','\r\n'));
fgetl(fid0); % Reading '#END SURFACE'
fclose(fid0);
%index = [1 2 3 4 5 6 7 8 11 9 10 12 13];
%ELEMENTS = ELEMENTS(:,index);


%WRITING FemurInput.mat
nyy = 0.3;
force = 1000;
nodes = NODES(:,1:4);


switch Method
    case 'A'
        %Write the Method A material file
        El_E = ELEMENTS(:,12);
    case 'B'
        %Write the Method B material file
        El_E = ELEMENTS(:,13);
end
elements = [ELEMENTS(:,1:5) El_E nyy*ones(length(El_E(:,1)),1)];
el_col = elements(:,2:5);
el_col = el_col(:);
I = intersect(nodes(:,1),el_col);
nodes = nodes(I,:);


%shaft support
MinShaftZ = min(nodes(:,4));
I = nodes(:,4) < MinShaftZ+5;
bc_nodes1 = nodes(I,1);
nn = length(bc_nodes1);
bc1 = [bc_nodes1 NaN(nn,1) NaN(nn,1) zeros(nn,1)];
    
%head support
MinHeadX = min(nodes(:,2));
I = nodes(:,2) < MinHeadX+5;
bc_nodes2 = nodes(I,1);
nn = length(bc_nodes2);
bc2 = [bc_nodes2 zeros(nn,1) NaN(nn,1) NaN(nn,1)];
    
%greater trochanter support
MaxGTX = max(nodes(:,2));
I = nodes(:,2) > MaxGTX-5;
bc_nodes3 = nodes(I,1);
nn = length(bc_nodes3);
bc3 = [bc_nodes3 -1*ones(nn,1) zeros(nn,1) NaN(nn,1)];

bcs = [bc1; bc2; bc3];

load1 = [1 0 0 0];
load2 = [];
save('FemurInput.mat','nodes','elements','bcs','load1','load2') 
    