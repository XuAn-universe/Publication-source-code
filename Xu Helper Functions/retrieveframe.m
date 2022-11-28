function im = retrieveframe(app, Exp_Path, trialID, timestamp)
readerobj1 = VideoReader([Exp_Path '\PG1camera' num2str(trialID) '.avi']);
readerobj2 = VideoReader([Exp_Path '\PG2camera' num2str(trialID) '.avi']);
video3location = strsplit(Exp_Path, filesep);
video3location{end-2} = [video3location{end-2} ' Camera3'];
video3location{end-1} = [video3location{end-1} ' PG3'];
video3location = strjoin(video3location, filesep);
readerobj3 = VideoReader([video3location '\PG3camera' num2str(trialID) '.avi']);

cam1track = find_tracking_filename(Exp_Path, ['PG1camera' num2str(trialID)]);
cam2track = find_tracking_filename(Exp_Path, ['PG2camera' num2str(trialID)]);
cam3track = find_tracking_filename(video3location, ['PG3camera' num2str(trialID)]);
filename = cell(1, 3);
filename{1} = cam1track;
filename{2} = cam2track;
filename{3} = cam3track;
num_track = cell(1, 3);
parfor i = 1:3 %
    if ~isempty(filename{i})
        [num_track{i}, ~, ~] = xlsread(filename{i});
    end
end
num_track1 = num_track{1};
num_track2 = num_track{2};
num_track3 = num_track{3};

currenttime = timestamp;
readerobj1.CurrentTime = currenttime-1/readerobj1.FrameRate;
im(:, :, :, 1) = get_frame_labels(app, readerobj1, num_track1);

readerobj2.CurrentTime = currenttime-1/readerobj2.FrameRate;
im(:, :, :, 2) = get_frame_labels(app, readerobj2, num_track2);

readerobj3.CurrentTime = currenttime-1/readerobj3.FrameRate;
im(:, :, :, 3) = get_frame_labels(app, readerobj3, num_track3);
end