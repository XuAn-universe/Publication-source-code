% generate color gradient as colorbar
%% for single color
nrow = 640;
cmap = ones(nrow, 1, 3);
cmap(:, :, 1) = 0.37; % or 0.68
cmap(:, :, 2) = (linspace(0, 1, nrow))';
im = imresize(hsv2rgb(cmap), [nrow 50], 'nearest');
figure;
imshow(flipud(im), []);

%% for gray scale
nrow = 640;
cmap = ones(nrow, 1, 3);
cmap(:, :, 2) = 0;
cmap(:, :, 3) = (linspace(0, 1, nrow))';
im = imresize(hsv2rgb(cmap), [nrow 50], 'nearest');
figure;
imshow(im, []);