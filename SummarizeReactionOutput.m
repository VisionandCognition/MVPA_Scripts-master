% DECODERRESULTS SVMWiedeField_StimDecoder

miceopt =  {'Alladin','Chief','Esmeralda','Frey'} %Happy'}%%% %options for mice
ReactionOpt = {'Error','Hit'};%{'Error','Hit','Miss'};
basel = [-300,-50];
TW = {[-300,-50],[100,500],[750,1350],[1350,1950]};
TWNames = {'Baseline','Visual','Delay','Response'};
nboot = 1000; 
spotanalysis = 0; %1 if spotlight, 0 if wide
DataFolder = 'I:\SARA\ReactionResults\';
OriDataPath =  ['\\vcnin\mouse_working_memory\Data4Class\TMPData\']; % Set path

posmap = fliplr([linspace(1,1,128);linspace(0,1,128);zeros(1,128)]);
% blackmap = fliplr([linspace(0.2,0.40,12);linspace(0.2,0.40,12);linspace(0.2,0.40,12)]);
negmap = fliplr([zeros(1,128);linspace(1,1,128);fliplr(linspace(0,1,128))]);
PSCOREMAP = fliplr(cat(2,posmap,negmap))';
blackval = round(0.95*size(PSCOREMAP,1)/2);
blackrange = (size(PSCOREMAP,1)/2)-blackval:(size(PSCOREMAP,1)/2)+(blackval-1);
blackmap = [fliplr(linspace(0,0.6,blackval)),linspace(0,0.6,blackval)]; %make 0.6 or sth instead of 1 to have more 'abrupt' black to color
PSCOREMAP(blackrange,:) = PSCOREMAP(blackrange,:).*repmat(blackmap,[3,1])';
for i = 1:3
    PSCOREMAP(:,i) = smooth(PSCOREMAP(:,i),5);
end
x = 1:256;
y = 1:256;
X = meshgrid(x,y);


ActSupColorMap = fliplr(cat(2,posmap,negmap))';



%% LOAD IN ORIGINAL DATA (non boot-strapped)
if spotanalysis
    AllBetaMaps = nan(length(miceopt),length(TW),262,262);

else
    AllBetaMaps = nan(length(miceopt),length(TW),800,800);
end
PERFORMANCE = nan(length(miceopt),length(TW));

 for midx = 1:length(miceopt)
        for twid = 1:length(TW)
            try
                tmp = load(fullfile(DataFolder,[miceopt{midx} '_time' num2str(TW{twid}(1)) '-' num2str(TW{twid}(2)) '.mat']));
                if spotanalysis                    
                    tmp.SVMMAP{1}(tmp.SVMMAP{1}==0 ) = nan;
                    AllBetaMaps(midx,twid,:,:) = tmp.SVMMAP{1};
                    PERFORMANCE(midx,twid) = nanmean(tmp.SVMMAP{1}(:));
                else
                    AllBetaMaps(midx,twid,:,:) = tmp.PERF{1}.beta;
                    PERFORMANCE(midx,twid) = nanmean(tmp.PERF{1}.svmperf);
                end
                
            catch ME
                try
                    if ~isempty(tmp.SVMMAP) && midx == 1 && ridx ==1 && twid == 1
                        AllBetaMaps = nan(length(miceopt),length(ReactionOpt),length(TW),size(tmp.SVMMAP{1},1),size(tmp.SVMMAP{1},2));
                    end
                    AllBetaMaps(midx,ridx,twid,:,:) = tmp.SVMMAP{1};
                catch ME
                    keyboard
                    disp(ME)
                    continue
                end
            end
        end
    BrainModel{midx} = load(fullfile(OriDataPath,miceopt{midx}, 'brainareamodel.mat'))
end
clear tmp

%Remove areas
throwawayareas = find(cellfun(@isempty,BrainModel{midx}.Model.Rnames));
throwawayareas = [throwawayareas; find(cellfun(@(X) ismember(X,{'OlfactoryBulb','fibrtracts','InfCol','SupColSens','ECT','Apo','Av'}),BrainModel{midx}.Model.Rnames))];

AreasOfInterest = BrainModel{midx}.Model.Rnames;
AreasOfInterest(throwawayareas) = [];

%% Figures & Summaries
BetaAvPerArea = nan(length(miceopt),length(TW),length(AreasOfInterest));
BetaStdPerArea = BetaAvPerArea;
for midx = 1:length(miceopt)
    tmp = abs(AllBetaMaps(midx,:,:,:));
    quantval1 = quantile(tmp(:),.05);
    quantval2 = quantile(tmp(:),.95);
    BETAFIG = figure('name',['BetaMaps' miceopt{midx}]);
    for twid = 1:length(TW)
        hh = subplot(1,length(TW),twid)
        if spotanalysis
            h = imagesc((squeeze(abs(AllBetaMaps(midx,twid,:,:)))),[0.1 1]);
            axis square
            set(h,'AlphaData',~isnan(squeeze(AllBetaMaps(midx,twid,:,:))))
             hold on
             plot(round(BrainModel{midx}.Model.AllX./(800/size(AllBetaMaps,3))), round(BrainModel{midx}.Model.AllY./(800/size(AllBetaMaps,4))),'k.')
        else
            h = imagesc(smooth2a(squeeze(abs(AllBetaMaps(midx,twid,:,:))),3),[quantval1 quantval2]);
            axis square
            set(h,'AlphaData',~isnan(squeeze(AllBetaMaps(midx,twid,:,:))))
             hold on
             plot(BrainModel{midx}.Model.AllX, BrainModel{midx}.Model.AllY,'k.')
            
        end
       
        title([TWNames{twid}])
        colormap(posmap')
    end
    
    Location = get(hh,'Position');
    colorbar
    set(hh,'Position',Location)

    saveas(BETAFIG,fullfile(DataFolder,['BetaMaps_' miceopt{midx} '.bmp']))
    saveas(BETAFIG,fullfile(DataFolder,['BetaMaps_' miceopt{midx} '.fig']))
    hgexport(BETAFIG,fullfile(DataFolder,['BetaMaps_' miceopt{midx} '.eps']))
    
    for areaid = 1:length(AreasOfInterest)
        %Extract brain region
        Borders = BrainModel{midx}.Model.Boundaries{strcmp(BrainModel{midx}.Model.Rnames,AreasOfInterest(areaid))};
        mask = false(size(AllBetaMaps,3),size(AllBetaMaps,4));
        for roi2dx = 1:length(Borders)
            mask2 = poly2mask(round(Borders{roi2dx}(:,1)./(800/size(mask,1))),round(Borders{roi2dx}(:,2)./(800/size(mask,2))),size(AllBetaMaps,3),size(AllBetaMaps,4));
            mask2 = bwmorph(mask2,'shrink',1);
            mask(mask2) = 1;
        end
        for ridx = 1:length(ReactionOpt)
            for twid = 1:length(TW)
                tmp = squeeze(AllBetaMaps(midx,twid,:,:));
                tmp(~mask) = nan;
                BetaAvPerArea(midx,twid,areaid) = nanmean(abs(tmp(:)));
                BetaStdPerArea(midx,twid,areaid) = nanstd(abs(tmp(:))./sum(~isnan(tmp(:))));
            end
        end
    end
    BARPLOT = figure('name',['BarPlots_' miceopt{midx}]);
    figure(BARPLOT)
    barwitherr(squeeze(BetaStdPerArea(midx,:,:))',squeeze(BetaAvPerArea(midx,:,:))')
    title(['Average Beta Values'])
    set(gca,'XTick',[1:length(AreasOfInterest)],'XTickLabel',AreasOfInterest)
    rotateXLabels(gca,45)
       
    legend(TWNames)
    
    saveas(BARPLOT,fullfile(DataFolder,['BarPLotsBetaMaps_' miceopt{midx} '.bmp']))
    saveas(BARPLOT,fullfile(DataFolder,['BarPLotsBetaMaps_' miceopt{midx} '.fig']))
    hgexport(BARPLOT,fullfile(DataFolder,['BarPLotsBetaMaps_' miceopt{midx} '.eps']))
    
    
    
end
%% Barplot
 AVBARS = figure('Name','AllAveragedBarPlots');
for twid = 1:length(TW)
    barwitherr(squeeze(nanstd(BetaAvPerArea(:,:,:),[],1))',squeeze(nanmean(BetaAvPerArea(:,:,:),1))')
    title(TWNames{twid})
    ylabel(['Average +/ std Beta Values'])
    set(gca,'XTick',[1:length(AreasOfInterest)],'XTickLabel',AreasOfInterest)
    rotateXLabels(gca,45)    
end
legend(TWNames)
saveas(AVBARS,fullfile(DataFolder,['BetaMapsAveragesBars.bmp']))
saveas(AVBARS,fullfile(DataFolder,['BetaMapsAveragesBars.fig']))
hgexport(AVBARS,fullfile(DataFolder,['BetaMapsAveragesBars.eps']))
   

