function Obj=SOFAconvertKosice2SOFA(root)
%SOFAconvertKosice2SOFA - converts from irs (Kosice format) to SOFA format
%   OBJ=SOFAconvertKosice2SOFA(root)
% 
%   SOFAconvertKosice2SOFA(root) converts BRIRs provided by Kosice stored in directory root to a SOFA object.
%
%   Input parameters:
%     root : directory containing the BRIRs (called nkce)
% 
%   Output parameters:
%     Obj : New SOFA object (SOFA format)

% #Author: Petra Kovacs (based on SOFAcovertTUBerlin2SOFA)
%
% Copyright (C) Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.2 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: https://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 


%% Get an empy conventions structure
Obj = SOFAgetConventions('SingleRoomSRIR');

%% Create empty matrix
M = 1;
Obj.SourcePosition = zeros(M,3);
Obj.Data.IR = zeros(M,2,32767);

%% Define elevation and sampling rate
ele = 0;
Obj.Data.SamplingRate = 44100; % based on Shinn-Cunningham et al., 2005, JASA

%% Read in and fill the data
% nkcePath = uigetdir;
% addpath(nkcePath);
% directory = dir(nkcePath);
addpath(root);
directory = dir(root);

for file = 3:length(directory)
    filename = directory(file).name;
    if strcmp(filename(1:4), "imp_")
        load(filename, 'impl', 'impr');
        azi = str2double(filename(10:12));
        distance = str2double(filename(14:16))/100; % converted from cm to m
        Obj.Data.IR(file-2,1,:) = shiftdim(impl,-2); % impl is [N M], data.IR must be [M R N]
        Obj.Data.IR(file-2,2,:) = shiftdim(impr,-2);
        Obj.SourcePosition(file-2,:) = [azi ele distance]; % When the spherical coordinate system is required, the format is (azimuth elevation distance)
    end
end

%% Fill with attributes
Obj.GLOBAL_ListenerShortName = 'KEMAR';
Obj.GLOBAL_History='Converted from the Kosice format';
Obj.GLOBAL_License = 'Creative Commons Attribution-NonCommercial-ShareAlike 3.0';
Obj.GLOBAL_ApplicationName = 'Distance-dependent BRIR from Kosice';
Obj.GLOBAL_ApplicationVersion = '1.0';
Obj.GLOBAL_AuthorContact = 'norbert.kopco@upjs.sk';
Obj.GLOBAL_References = ['N. Kopco, K.K. Doreswamy, S. Huang, S. Rossi, J. Ahveninen. ', ...
    'Cortical auditory distance representation based on direct-to-reverberant',...
    'energy ratio. In Neuroimage, 208 (2020)'];
Obj.GLOBAL_DatabaseName = 'Kosice';
Obj.GLOBAL_Title = 'BRIR';
Obj.GLOBAL_Organization = 'P. J. Safarik University, Kosice, Slovakia';

%% Fill the mandatory variables
% Obj.ListenerPosition = [0 0 0]; % has to be the same size as source position (29*3)
Obj.ListenerPosition = zeros(size(Obj.SourcePosition));
Obj.ListenerView = [1 0 0];
Obj.ListenerUp = [0 0 1];
Obj.Data.Delay = [0 0];
  
%% Save
SOFAsave('KosiceBRIR.sofa',Obj);
