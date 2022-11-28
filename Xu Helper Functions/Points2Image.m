function Points2Image(app, FrameRate, readerobj, num_track, csvdata)
trial = app.TrialsListBox.Value;
if numel(trial) ~= 1
    helpdlg('You can only choose one trial');
    return;
end
if isempty(num_track)
    return;
end

value = app.Slider.Value;
time_range = [app.timerangesEditField.Value app.toEditField.Value];
label_table = table2array(app.UITable.Data);
labelID = [3 10 16];
color_NI = [0 0 0.5; 0.5 0 0; 0 0.5 0];
color_I = [0 0 1; 1 0 0; 0 1 0];
brightness = app.BrightnessEditField.Value;

if strcmp(readerobj.Name(end-2:end), 'avi')
    readerobj.CurrentTime = value-1/readerobj.FrameRate;
elseif strcmp(readerobj.Name(end-2:end), 'mp4')
    readerobj.CurrentTime = value;
end
im = readFrame(readerobj);
if brightness ~= 1
    im = rgb2hsv(im);
    intensity = im(:, :, 3);
    intensity = mat2gray(intensity, [0 brightness]);
    im(:, :, 3) = intensity;
    im = hsv2rgb(im);
end
figure;
imshow(im);
hold on;
t = (1:size(num_track, 1))'/FrameRate;
for i = 1:numel(labelID)
    x = num_track(:, labelID(i)*3-1);
    y = num_track(:, labelID(i)*3);
    Istart = t(find(csvdata(:, 1) == 1, 1));
    Istop = t(find(csvdata(:, 1) == 1, 1, 'last'));
    
    x_NI = x(t >= time_range(1) & t <= time_range(2) & t < Istart & num_track(:, labelID(i)*3+1) > label_table(labelID(i), 1));
    y_NI = y(t >= time_range(1) & t <= time_range(2) & t < Istart & num_track(:, labelID(i)*3+1) > label_table(labelID(i), 1));
    
    x_I = x(t >= time_range(1) & t <= time_range(2) & csvdata(:, 1) == 1 & num_track(:, labelID(i)*3+1) > label_table(labelID(i), 1));
    y_I = y(t >= time_range(1) & t <= time_range(2) & csvdata(:, 1) == 1 & num_track(:, labelID(i)*3+1) > label_table(labelID(i), 1));
    
    quiver(x_NI, y_NI, diff([x_NI; x_I(1)]), diff([y_NI; y_I(1)]), 0, 'Color', color_NI(i, :));
    quiver(x_I(1:end-1), y_I(1:end-1), diff(x_I), diff(y_I), 0, 'Color', color_I(i, :));
    
    x_NI = x(t >= time_range(1) & t <= time_range(2) & t > Istop & num_track(:, labelID(i)*3+1) > label_table(labelID(i), 1));
    y_NI = y(t >= time_range(1) & t <= time_range(2) & t > Istop & num_track(:, labelID(i)*3+1) > label_table(labelID(i), 1));
    
    quiver([x_I(end); x_NI(1:end-1)], [y_I(end); y_NI(1:end-1)], diff([x_I(end); x_NI]), diff([y_I(end); y_NI]), 0, 'Color', color_NI(i, :));
end
legend('Nose Before', 'Nose During', 'Nose After', 'Left Paw Before', 'Left Paw During', 'Left Paw After', 'Right Paw Before', 'Right Paw During', 'Right Paw After');