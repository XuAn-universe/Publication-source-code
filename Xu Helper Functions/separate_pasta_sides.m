function [IDleft_front, IDright_front, IDleft_top, IDright_top] = separate_pasta_sides(ctop, cbottom)
%*---------------------------------------------------------------------*
% Huang Lab
% Duke University
% Author : Xu An, Dec 2021
% xu.an@duke.edu
% Version: 1.0
%*---------------------------------------------------------------------*
IDleft_front = [];
IDright_front = [];
IDleft_top = [];
IDright_top = [];
nbite = size(ctop, 3);
ctop = squeeze(ctop);
cbottom = squeeze(cbottom);
if nbite == 1
    ctop = ctop';
    cbottom = cbottom';
end
for i = 1:size(ctop, 2)
    if ~any(isnan([ctop([1 3], i); cbottom([1 3], i)]))
        if ctop(3, i) >= cbottom(3, i)
%             if (ctop(1, i)-cbottom(1, i) <= 0 && cbottom(1, i) > -5) || ctop(1, i) > 5
            if (ctop(1, i)-cbottom(1, i) <= 0 && ctop(1, i) > 2)  || cbottom(1, i) > 7
                IDleft_front = [IDleft_front i];
                if ~any(isnan([ctop(2, i); cbottom(2, i)]))
                    IDleft_top = [IDleft_top i];
                end
            else
                IDright_front = [IDright_front i];
                if ~any(isnan([ctop(2, i); cbottom(2, i)]))
                    IDright_top = [IDright_top i];
                end
            end
        else
%             if (cbottom(1, i)-ctop(1, i) <= 0 && ctop(1, i) > -5) || cbottom(1, i) > 5
            if (cbottom(1, i)-ctop(1, i) <= 0 && cbottom(1, i) > 2)  || ctop(1, i) > 7
                IDleft_front = [IDleft_front i];
                if ~any(isnan([ctop(2, i); cbottom(2, i)]))
                    IDleft_top = [IDleft_top i];
                end
            else
                IDright_front = [IDright_front i];
                if ~any(isnan([ctop(2, i); cbottom(2, i)]))
                    IDright_top = [IDright_top i];
                end
            end
        end
    end
end