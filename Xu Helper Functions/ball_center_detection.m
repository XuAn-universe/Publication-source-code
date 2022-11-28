% Detect a same ball using multiple cameras to make the center of the ball as the common reference point
% By Xu An 11/9/2019 @Cold Spring Harbor Laboratory
% xan@cshl.edu
%%
figure;
imshow(Ball_Cam1);
d = imdistline;
%%
delete(d);
%%
[centers1, radii1] = imfindcircles(Ball_Cam1, [75 85], 'ObjectPolarity', 'bright', 'Sensitivity', 0.95);
h = viscircles(centers1, radii1);

%%
figure;
imshow(Ball_Cam2);
d = imdistline;
%%
delete(d);
%%
[centers2, radii2] = imfindcircles(Ball_Cam2, [85 95], 'ObjectPolarity', 'bright', 'Sensitivity', 0.95);
h = viscircles(centers2, radii2);

%%
figure;
imshow(Ball_Cam3);
d = imdistline;
%%
delete(d);
%%
[centers3, radii3] = imfindcircles(Ball_Cam3, [70 80], 'ObjectPolarity', 'bright', 'Sensitivity', 0.95);
h = viscircles(centers3, radii3);

%%
ref_point = [centers1; centers2; centers3];