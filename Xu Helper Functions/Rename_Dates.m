%%
DataFolder = 'D:\Fiber Photometry Data';
year_range = 15:21;
DirList = dir(DataFolder);
for n = 3:size(DirList)
    ExpName = DirList(n).name;
    if ismember(str2double(ExpName(5:6)), year_range)
        movefile([DataFolder '\' ExpName], [DataFolder '\20' ExpName(5:6) ExpName(1:4) ExpName(7:end)]);
    end
end
disp('Done !');