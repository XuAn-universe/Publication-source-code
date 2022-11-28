function Update_Labels4LickTrials(app, Exp_Path)
%*---------------------------------------------------------------------*
% Huang Lab
% Duke University
% Author : Xu An, Jan 2022
% xu.an@duke.edu
% Version: 1.0
%*---------------------------------------------------------------------*
value = app.TrialsListBox.Value;
value = sort(value);
trials = numel(value);

for i = 1:trials
    temp = load([Exp_Path '\LabelledEvents' num2str(value(i)) '.mat']);
    LabelledEvents = temp.LabelledEvents;
    if isempty(LabelledEvents.BiteBoutStart)
        MouthOpen = LabelledEvents.MouthOpen;
        MouthClosed = LabelledEvents.MouthClosed;
        if ~isempty(MouthOpen) && ~isempty(MouthClosed)
            LabelledEvents.RetrievalStart = max(MouthOpen);
            LabelledEvents.MouthRetrievalStart = max(MouthOpen);
            LabelledEvents.MouthRetrievalEnd = max(MouthClosed);
            save([Exp_Path '\LabelledEvents' num2str(value(i)) '.mat'], 'LabelledEvents');
        end
    end
end
msgbox('Done !');