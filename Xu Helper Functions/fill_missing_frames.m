%%
pathname = {'I:\Free-moving Feeding Data\20221010_AngelHair_Fezf2CreER;PlexinD1flp_Gcamp-R1\exp001',...
    'I:\Free-moving Feeding Data\20221010_AngelHair_Fezf2CreER;PlexinD1flp_Gcamp-R2\exp001',...
    'I:\Free-moving Feeding Data\20221013_AngelHair_Fezf2CreER;PlexinD1flp_Gcamp-R1\exp001',...
    'I:\Free-moving Feeding Data\20221013_AngelHair_Fezf2CreER;PlexinD1flp_Gcamp-R2\exp001',...
    'I:\Free-moving Feeding Data\20221014_AngelHair_Fezf2CreER;PlexinD1flp_Gcamp-R1\exp001',...
    'I:\Free-moving Feeding Data\20221014_AngelHair_Fezf2CreER;PlexinD1flp_Gcamp-R2\exp001',...
    'I:\Free-moving Feeding Data Camera3\20221010_AngelHair_Fezf2CreER;PlexinD1flp_Gcamp-R1 PG3\exp001',...
    'I:\Free-moving Feeding Data Camera3\20221010_AngelHair_Fezf2CreER;PlexinD1flp_Gcamp-R2 PG3\exp001',...
    'I:\Free-moving Feeding Data Camera3\20221013_AngelHair_Fezf2CreER;PlexinD1flp_Gcamp-R1 PG3\exp001',...
    'I:\Free-moving Feeding Data Camera3\20221013_AngelHair_Fezf2CreER;PlexinD1flp_Gcamp-R2 PG3\exp001',...
    'I:\Free-moving Feeding Data Camera3\20221014_AngelHair_Fezf2CreER;PlexinD1flp_Gcamp-R1 PG3\exp001',...
    'I:\Free-moving Feeding Data Camera3\20221014_AngelHair_Fezf2CreER;PlexinD1flp_Gcamp-R2 PG3\exp001'};

thr = 1;

curpwd = pwd;
workbar(0, 'Processing...', 'Progress');
for i = 1:numel(pathname)
    cd(pathname{i});
    DirList = dir;
    total = length(DirList);
    for n = 3:total
        if strcmp(DirList(n).name(end-3:end), '.avi')
            csvname = [DirList(n).name(1:3) 'gpio' DirList(n).name(10:end-4) '.csv'];
            csvdata = load(csvname);
            if max(diff(csvdata(:, 2))) > 1
                lastid = nan;
                framecount = 0;
                readerobj = VideoReader(DirList(n).name);
                FrameRate = readerobj.FrameRate;
                vidObj = VideoWriter(DirList(n).name(1:end-4), 'MPEG-4');
                set(vidObj, 'FrameRate', FrameRate, 'Quality', 100);
                open(vidObj);
                while hasFrame(readerobj)
                    vidFrame = readFrame(readerobj);
                    if thr ~= 1
                        vidFrame = rgb2hsv(vidFrame);
                        intensity = vidFrame(:, :, 3);
                        intensity = mat2gray(intensity, [0 thr]);
                        vidFrame(:, :, 3) = intensity;
                        vidFrame = hsv2rgb(vidFrame);
                    end
                    framecount = framecount+1;
                    currentid = csvdata(framecount, 2);
                    if ~isnan(lastid)
                        framesmissing = currentid-lastid-1;
                        if framesmissing ~= 0
                            for j = 1:framesmissing
                                writeVideo(vidObj, vidFrame.*0);
                            end
                        end
                    end
                    writeVideo(vidObj, vidFrame);
                    lastid = currentid;
                end
                close(vidObj);
                clear readerobj;
                delete(DirList(n).name);
                
                csvdata_copy = nan(csvdata(end, 2)-csvdata(1, 2)+1, 2);
                csvdata_copy(:, 2) = 1:size(csvdata_copy);
                csvdata_copy(csvdata(:, 2)-csvdata(1, 2)+1, 1) = csvdata(:, 1);
                csvdata_copy(:, 1) = fillmissing(csvdata_copy(:, 1), 'nearest');
                dlmwrite(csvname, csvdata_copy, 'precision', 6);
            end
            disp(DirList(n).name);
        end
    end
    workbar(i/numel(pathname), pathname{i}, 'Progress');
end
cd(curpwd);
msgbox('Done');