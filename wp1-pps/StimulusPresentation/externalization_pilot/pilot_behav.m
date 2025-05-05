function pilot_behav

% Convert externalization pilot data to csv
%
% Input mat file:
% ext_resp_SUBNUM.mat
% - results.participant - participant number
% - results.cols - column names (except particpant number)
% - results.res - data (NTrials X NVariables, i.e. 128 X 6)
%
% Output csv columns:
% participant, source intensity condition, f0 condition, f0, azimuth,
% distance, externalization rating
% Size: NSubjects*NTrials X NVariables, i.e. 8*128 X 7
% 

dirName = uigetdir;
dirFile = dir(dirName);
numFile = size(dirFile);
addpath(dirName);

for i=3:numFile(1)
    filename = dirFile(i).name;
    load(filename, "results");
    
    sourceInt = results.res(:,1);
    f0cond = results.res(:,2);
    f0 = results.res(:,3);
    azimuth = results.res(:,4);
    distance = results.res(:,5);
    extRating = results.res(:,6);
    subNum = repmat(results.participant,length(sourceInt),1);
    
    varNames = {'subNum','sourceInt','f0cond','f0','azimuth','distance','extRating'};
    M = table(subNum, sourceInt, f0cond, f0, azimuth, distance, extRating, 'VariableNames',varNames);
    
    if i < 4
        writetable(M, 'externData.csv',"WriteMode","append","WriteVariableNames",true);
    else
        writetable(M, 'externData.csv',"WriteMode","append","WriteVariableNames",false);
    end

end