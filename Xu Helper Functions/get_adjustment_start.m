function tstart = get_adjustment_start(rstart, rend, lstart, lend)
tstart = [];
overlapID_all = [];
for i = 1:numel(rstart)
    overlapID = find((lstart <= rstart(i) & lend >= rstart(i)) | (lstart > rstart(i) & lstart <= rend(i)));
    tstart = [tstart; min([rstart(i); lstart(overlapID)])];
    overlapID_all = [overlapID_all; overlapID];
end
lstart(overlapID_all) = [];
tstart = [tstart; lstart];