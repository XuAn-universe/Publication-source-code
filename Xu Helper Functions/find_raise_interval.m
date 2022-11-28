function result = find_raise_interval(statepath)
%*---------------------------------------------------------------------*
% Huang Lab
% Duke University
% Author : Xu An, Jan 2022
% xu.an@duke.edu
% Version: 1.0
%*---------------------------------------------------------------------*
chewID = []; % start of chew
handleID = []; % end of handle
chew2handleID = []; % start of handle
chewtemp = [];
handletemp = [];
for i = 1:numel(statepath)
    if i == 1
        if statepath(i) == 1
            chewtemp = [chewtemp i];
        else
            handletemp = [handletemp i];
        end
        laststate = statepath(i);
    else
        if statepath(i) == laststate
            if statepath(i) == 1
                chewtemp = [chewtemp i];
                if i == numel(statepath)
                    chewID = [chewID min(chewtemp)];
                end
            else
                handletemp = [handletemp i];
                if i == numel(statepath)
                    handleID = [handleID max(handletemp)];
                end
            end
        else
            if laststate == 1
                chewID = [chewID min(chewtemp)];
                chewtemp = [];
                handletemp = [handletemp i];
                if i == numel(statepath)
                    handleID = [handleID max(handletemp)];
                end
                chew2handleID = [chew2handleID i];
            else
                handleID = [handleID max(handletemp)];
                handletemp = [];
                chewtemp = [chewtemp i];
                if i == numel(statepath)
                    chewID = [chewID min(chewtemp)];
                end
            end
            laststate = statepath(i);
        end
    end 
end
if ~isempty(chewID) && ~isempty(handleID)
    if chewID(1) > handleID(1)
        handleID(1) = [];
    end
    if numel(chewID) > numel(handleID)
        chewID(end) = [];
    end
end
if ~isempty(chewID) && ~isempty(handleID)
    result = [chewID; handleID; chew2handleID];
    result = floor(result);
else
    result = [];
end