 MouseOpt = {'Innoko'};
 timeopt = [1:4];
%Run SVM on Desktop
for Mouseidx = 1:length(MouseOpt)
    parfor TWidx = 1:length(timeopt)        
        SVMWideField_ReactionDecoder(Mouseidx,TWidx)
    end
end
