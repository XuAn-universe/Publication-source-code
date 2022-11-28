function im = get_frame_labels(app, readerobj, num_track)
im = readFrame(readerobj);
dotsize = 3;
alphavalue = 0.7;
labelcolor = jet(256)*2^8;
frameID = round(readerobj.CurrentTime*readerobj.FrameRate);
label_table = table2array(app.UITable.Data);
nlabels = size(label_table, 1);
for i = 1:nlabels
    if label_table(i, 2) && num_track(frameID, i*3+1) > label_table(i, 1)
        im = insertShape(im, 'FilledCircle', [num_track(frameID, i*3-1) num_track(frameID, i*3) dotsize], 'Color',...
            labelcolor(round((i-1)*255/(nlabels-1))+1, :), 'Opacity', alphavalue);
    end
end
end