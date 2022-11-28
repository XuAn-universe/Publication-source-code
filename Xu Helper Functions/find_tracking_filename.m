function filename = find_tracking_filename(Exp_Path, camID)
filename = [];
length_camID = length(camID);
curpwd = pwd;
cd(Exp_Path);
DirList = dir;
for i = 3:size(DirList)
    if numel(DirList(i).name) > length_camID
        if strcmp(DirList(i).name(1:length_camID), camID) && strcmp(DirList(i).name(end-3:end), '.csv') &&...
                ~isempty(strfind(DirList(i).name, 'DeepCut'))
            filename = [Exp_Path '\' DirList(i).name];
        end
    end
end
cd(curpwd);
end