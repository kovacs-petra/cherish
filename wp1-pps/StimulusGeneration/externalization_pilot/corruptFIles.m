corrupt = stimStruct.stim(:,6,1);
sizes=size(corrupt);
for ii = 1:sizes(1)
    % for jj = 1:sizes(2)
        for hh = 1:sizes(2)
            cellG = corrupt{ii,hh};
            figure; plot(cellG{1});
            tit = strcat('f0=',num2str(ii), ' azi=',num2str(hh));
            title(tit);
        end
    % end
end