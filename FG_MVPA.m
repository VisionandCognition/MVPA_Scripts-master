function FG_MVPA(Mouseidx,TWidx,roidx,roimethod)
% Decoding --> Use MVPA (SVM) to decode orientation of the stimulus for one
% session
if ~exist('Mouseidx','var')
    Mouseidx =1;
    TWidx = 1;
    roidx = 1;
    roimethod = 3;
end

if ischar(Mouseidx)
    Mouseidx = str2num(Mouseidx);
end
if ischar(TWidx)
    TWidx = str2num(TWidx);
end
if ischar(roidx)
    roidx = str2num(roidx);
end
if ischar(roimethod)
    roimethod = str2num(roimethod);
end

%Station
clear notused
%% User Defined Input
miceopt = {'Jules','Marsellus','Vincent','Marsellus'} %options for mice
Stim2Check = 'FGTask'%Name of the stimulus as written in the LOG-file
%Timelimit: Don't need data from time after this.
timelimit1 = 2500; %ms
BGOpt = [0,1];
BGOptNames = {'Grey','Contrast'};
if strcmp(Stim2Check,'FGTask')
    TW = {[-300,0],[100,500]};
    TWNames = {'Baseline','Visual'};
end
cd
if ispc
    DataPath = '\\vcnin\mouse_working_memory\Data4Class\TMPData'; % Set path
    ResultPath = 'I:\SARA\TMPResults'; % Set path
    ScriptsPath = 'I:\SARA\MVPA_Scripts\FG'; %Set path
    TWidx = 2;
else
    DataPath = fullfile(cd,miceopt{Mouseidx}) % Set path
    ResultPath = fullfile(cd,'TMPResults') % Set path
    ScriptsPath = fullfile(cd,'MVPA_Scripts') %Set path
end

if ~exist(ResultPath,'dir')
    mkdir(ResultPath)
end

addpath(genpath(DataPath))
addpath(genpath(ResultPath))
% addpath(genpath(ScriptsPath))
nfolds = 10; %Number of folds for crossvalidation & nested cross-validation
n = -1:0.50:9; %Range of what parameter c could be (c = 2^n(i))
PredictorTimerange = [-100 2500]; %Timerange lower and upperbound for which predictions are done (in ms)

tempstation = fullfile(ResultPath,'TMPMatlab');
if ~exist(tempstation,'dir')
    mkdir(tempstation)
end
%Load info
load(fullfile(DataPath,'sessionstruct.mat'))
load(fullfile(DataPath,'randvec.mat'))

if ispc    
    DataPath = ['\\vcnin\mouse_working_memory\Data4Class\TMPData\' miceopt{Mouseidx}]; % Set path
end

%% Reading datapaths etc.
paths = info.paths;
logs = info.logs;

mousecount = 0;
for midx = Mouseidx %For this mouse
    if sum(~cellfun(@isempty, {logs{midx,:,:}})) < 1 %If not recorded that day, skip
        continue
    end
    mousecount = mousecount+1;
    mouse = miceopt{midx};
    sessioncount = 0;
    clear RawData
    
    %% Take models
    %     Alan Brain
    if roimethod == 1
        load(fullfile(DataPath,'brainareamodel.mat'))
        regio2take = find(~strcmp(Model.Rnames,''));
    elseif roimethod == 2
        %     OR
        % ROIs
        load(fullfile(DataPath,[mouse 'EvokedActivROIs']))
        %make a model of this
        regio2take = [1:length(rois)]';
        Model.Boundaries = cell(1,length(rois));
        for ii = 1:length(rois)
            Model.Boundaries{ii}{1} = [rois{ii}.xi rois{ii}.yi];
            Model.Rnames{ii} = ['ROI' num2str(ii)];
        end
    elseif roimethod == 3
        load(fullfile(DataPath,'brainareamodel.mat'))
        regio2take = 1;        
    end
        
    if roidx > length(regio2take)
       disp(['too large number for the number of rois that exist'])
       return
    end
    if roidx == length(regio2take)
        disp(['Check whether all rois are included..'])
    end

    
    delete(fullfile(tempstation,[mouse '_' Model.Rnames{regio2take(roidx)} '_TW' num2str(TWidx) '_tmpfile.mat']))
    TMPMAT = matfile(fullfile(tempstation,[mouse '_' Model.Rnames{regio2take(roidx)} '_TW' num2str(TWidx) '_tmpfile.mat'])); %Make a workable matfile
    
    for didx = 1:size(logs,2) %Loop over days
        if sum(~cellfun(@isempty, {logs{midx,didx,:}})) < 1 %If not recorded that day, skip
            continue
        end
        for sidx = 1:size(logs,3) %If no xth session, continue
            if sum(~cellfun(@isempty,{logs{midx,didx,sidx}})<1)
                continue
            end
            sessioncount = sessioncount+1;
            clear tosave;
            clear LOG
            clear this
            tmppath = paths{midx,didx,sidx};
            date = strsplit(tmppath,mouse);
            date = date{3}(1:end-1) %Find date
            expnr = strsplit(tmppath,mouse);
            expnr = str2num(expnr{end});%find session nr
            
            disp(['Loading data ' mouse ', day ' date ', session ' num2str(expnr)])
            
            %% Log file
            if ispc
                load(fullfile('\\vc2nin\WBImaging\',mouse,[mouse date],[mouse num2str(expnr)],[mouse num2str(expnr) '.mat']));
            else
                load(fullfile(DataPath,[mouse date],[mouse num2str(expnr)],[mouse num2str(expnr) '.mat']));
            end
            
            if exist('tosave','var')
                try
                    LOG=tosave.LOG;
                catch
                    LOG=tosave.Log;
                end
            end
            if exist('Log','var')
                LOG = Log;
                clear Log;
            end
            
            if strcmp(Stim2Check,'DelayedOriTuningSound')
                %Make this.log.Orientation longer with nans
                LOG.Orientation(end:length(LOG.Reaction)) = 500; %Only goes till 360
                
                while ~isfield(LOG,'correctReaction') %Check whether reactions were registered okay
                    CheckReactions(fullfile('\\vc2nin\WBImaging\',mouse,[mouse date],[mouse num2str(expnr)],[mouse num2str(expnr) '.mat']))
                    tmp = load([folder expname '\' mouse expnum '.mat']);
                    if isfield(tmp,'tosave')
                        tmp = tmp.tosave;
                    end
                end
                LOG.Reaction = LOG.correctReaction; %Change the reactions into checked reactions
            end
            %% Create timewindows
            OriOpt = unique(LOG.Orientation);
            if isfield(LOG,'Side')
                SideOpt = unique(LOG.Side);
            else
                SideOpt = 1;
            end
            if ~iscell(SideOpt)
                SideOpt = {num2str(SideOpt)};
            end
            if isfield(LOG,'Reactions') || isfield(LOG,'Reaction')
                ReactionOptTMP = {'Miss','Hit','Error','Too Early','TooFast'};
                LOG.Condition = zeros(length(LOG.Reaction), 1);
            end
            
            count = 0;
            for oidx = 1:length(OriOpt)
                for soidx = 1:length(SideOpt)
                    if isfield(LOG,'Reactions') | isfield(LOG,'Reaction') %active
                        for rtmpidx = 1:length(ReactionOptTMP)
                            count = count + 1;
                            LOG.Condition(strcmp(LOG.Reaction,ReactionOptTMP{rtmpidx})& LOG.Orientation == OriOpt(oidx) & ...
                                strcmp(LOG.Side,SideOpt{soidx})) = count;
                            ConditionNames{count} = [ReactionOptTMP{rtmpidx} ' Ori' num2str(OriOpt(oidx)) ' Side ' SideOpt{soidx}];
                        end
                    else %Passive
                        count = count + 1;
                        LOG.Condition(LOG.Orientation == OriOpt(oidx) & ...
                            strcmp(LOG.Side,SideOpt{soidx})) = count;
                        ConditionNames{count} = ['Ori' num2str(OriOpt(oidx)) ' Side ' SideOpt{soidx}];
                    end
                    
                end
            end
            
            LOG.Conditions = unique(LOG.Condition);
            cvec = LOG.Conditions;
            if size(cvec,1) > size(cvec,2)
                cvec = cvec'
            end
            cvec(cvec==0) = [];
            ConditionNames(cvec)
            
            idx = find(~cellfun(@isempty,strfind(ConditionNames,'Too Early')));
            ConditionNames(idx) = cellfun(@(X) strrep(X,X(strfind(X,'Too Early'):9),'TooEarly'),ConditionNames(idx),'UniformOutput',0);
            
            %Average over orientations
            conditionparts = cellfun(@(X) strsplit(X,' '),ConditionNames(cvec),'UniformOutput',0);
            
            %Find all reactions
            reaction = cellfun(@(X) X{1},conditionparts,'UniformOutput',0); %Reaction
            orientation = cellfun(@(X) X{2},conditionparts,'UniformOutput',0); %orientations
            OriOpt = unique(orientation);
            side = cellfun(@(X) X{4},conditionparts,'UniformOutput',0); %SIdes
            SideOpt = unique(side);
            
        
            
            %%  Load 'drift correction'
            load(fullfile(DataPath,[mouse date],[mouse num2str(expnr)],'BASELINEMAT.mat'))
            
            %% Load data
            rawdatfiles = dir(fullfile(DataPath,[mouse date],[mouse num2str(expnr)],[mouse num2str(expnr) '_RawData*']));
            
            %Load data movement matrix
            load(fullfile(DataPath,[mouse date],[mouse num2str(expnr)],'ThrowAwayIdx.mat'))
            for ridx = 1:length(BGOpt)
                
                fulldelay = find(LOG.BGContrast==BGOpt(ridx) & LOG.Gavepassive == 0 & LOG.Ezbox == 0);
                for stidx = 1:length(SideOpt)
                    ccidx = find(strcmp(side,SideOpt{stidx}));
                    if isempty(ccidx)
                        nrtotal{sessioncount,stidx,ridx} = nrtotal{sessioncount-1,stidx,ridx};
                        continue
                    end
                    cidxcount = 0;

                    for cidx = ccidx
                        disp(['Loading data condition ' num2str(cidx) ' of ' num2str(length(cvec))])
                        clear conddata
                        cidxcount = cidxcount+1;
                        try
                            load(fullfile(DataPath,[mouse date],[mouse num2str(expnr)],rawdatfiles(strcmp({rawdatfiles(:).name},[mouse num2str(expnr) '_RawData_C' num2str(cidx) '.mat'])).name));
                            %                         keepidx = ismember(ctrials{cidx},fulldelay);
                            keepidx = ~removeidx(1:size(conddata,4),cidx)'&ismember(ctrials{cidx},fulldelay);
                            %Throw out the motion trials
                            if sessioncount==1 && cidxcount == 1
                                nrtotal{sessioncount,stidx,ridx} = sum(keepidx);
                            elseif cidxcount ~= 1
                                nrtotal{sessioncount,stidx,ridx} = nrtotal{sessioncount,stidx,ridx}+sum(keepidx);
                            else
                                nrtotal{sessioncount,stidx,ridx} = nrtotal{sessioncount-1,stidx,ridx}+sum(keepidx);
                            end
                            
                            if sum(keepidx)==0
                                continue
                            end
                            
                            m = whos(TMPMAT);
                            m = {m(:).name};
                            
                            tmp = single(conddata(:,:,(timeline>=TW{TWidx}(1) & timeline <= TW{TWidx}(2)),keepidx));
                            tmp(tmp==0)=nan;
                            tmp = squeeze(nanmean(tmp,3));
                            
                            %Cut out brain region
%                             Borders = Model.Boundaries{regio2take(roidx)};
%                             mask = zeros(size(tmp,1),size(tmp,2));
%                             for roi2dx = 1:length(Borders)
%                                 masktmp = poly2mask(Borders{roi2dx}(:,1),Borders{roi2dx}(:,2),size(tmp,1),size(tmp,2));
%                                 %Shrink to not have border effects
%                                 masktmp = bwmorph(masktmp,'shrink',1);
%                                 mask(masktmp==1)=1;
%                             end
%                             tmp(~repmat(mask,[1,1,size(tmp,3)])) = nan;
                            tmpname = [BGOptNames{ridx} '_' SideOpt{stidx}];
                            if ~ismember(tmpname,m)
                                if length(size(tmp))<3
                                    tmp = cat(3,tmp,nan(size(tmp)));
                                end
                                for i = 1:100:size(conddata,1)
                                    RawTmp(i:i+99,:,:) =  tmp(i:i+99,:,:)./BASELINEMAT(i:i+99,:,ctrials{cidx}(keepidx));          %baselinedrift
                                    %                             RawData{sessioncount,stidx}(i:i+99,:,:,:) =  tmp(i:i+99,:,:,:)./ permute(repmat(BASELINEMAT(i:i+99,:,ctrials{cidx}(keepidx)),[1,1,1,size(tmp,3)]),[1,2,4,3]);          %baselinedrift
                                end
                                eval(['TMPMAT.' tmpname ' = RawTmp;']);
                                clear RawTmp
                            elseif cidxcount ~= 1
                                 for i = 1:100:size(conddata,1)
                                    eval(['TMPMAT.' tmpname '(i:i+99,:,nrtotal{sessioncount,stidx,ridx}-sum(keepidx)+1:nrtotal{sessioncount,stidx,ridx}) = tmp(i:i+99,:,:)./BASELINEMAT(i:i+99,:,ctrials{cidx}(keepidx));'])
                                 end                                
                            else
                                for i = 1:100:size(conddata,1)
                                    eval(['TMPMAT.' tmpname '(i:i+99,:,nrtotal{sessioncount-1,stidx,ridx}+1:nrtotal{sessioncount,stidx,ridx}) = tmp(i:i+99,:,:)./BASELINEMAT(i:i+99,:,ctrials{cidx}(keepidx));'])
                                end
                            end
                        catch ME
                            disp(ME)
                            for i = 1:length(ME)
                                disp(ME.stack(i))
                            end
                            keyboard
                        end
                    end
                    
                end
            end
            timeline(timeline>timelimit1) = [];
            
            
        end
       
    end
    clear conddata
    clear SErrorval
    clear RawTmp
    
    nrtrials = nrtotal(sessioncount,:,:);
    nrtrials = [nrtrials{:}];
    %% MVPA - check for every brain region whether it's pixels can predict orientation of the stimulus.
    takenrt = min(nrtrials); %nr of trials to take (should be even number for both sides)
    
    if takenrt < 5
        disp(['Less than 10 trials to take... Skipping this mouse'])
        continue
    end
    %Make dataset for this reaction
    clear orilabelsnr
    clear TmpDat
    clear tmp
    clear BASELINEMAT
    clear randidx
    
    %Random vector of trials from trainingset for inner loop
    randomtidx = shamtrialvec(find(shamtrialvec<=length(nrtrials)*takenrt));      
    chunksz = floor(length(randomtidx)./nfolds); %Even chunks of data as nested-cross
    
    %% Start Decoding
    % Add paths necessary
    %Call TMPMAT
    TMPMAT = matfile(fullfile(tempstation,[mouse '_' Model.Rnames{regio2take(roidx)}  '_TW' num2str(TWidx) '_tmpfile.mat'])); %Make a workable matfile
    
    % Write data into GP understandable 'structs'
    m = whos(TMPMAT);
    m = {m(:).name};
    
    xpix = size(eval(['TMPMAT.' m{1}]),1);
    ypix = size(eval(['TMPMAT.' m{1}]),2);
    
    %Find non-nan ranges
    takethese = ~isnan(nanmean(eval(['TMPMAT.' m{1}]),3));
    takethese = find(takethese);
    nrpixels2take = length(takethese(:));
    
    tmp = nanmean(eval(['TMPMAT.' m{1}]),3);
    [r c] = find(~isnan(tmp));
    
    ReactionNumbers = 1:length(unique(BGOpt));
    SideNumbers = 1:length(unique(SideOpt));
    
    XDat = zeros(length(m)*takenrt,nrpixels2take,'single');
    YDat = zeros(2,length(m)*takenrt,'single');
    for mid = 1:length(m)
        tmp =  eval(['TMPMAT.' m{mid}]);
        randomtake = shamtrialvec(find(shamtrialvec<=size(tmp,3))  );
        tmp = tmp(:,:,randomtake(1:takenrt));
        for trialidx = 1:takenrt
            tmp2 = reshape(tmp(:,:,trialidx),[size(tmp,1)*size(tmp,2),1]);
            XDat((mid-1)*takenrt+trialidx,:) = tmp2(takethese,:);
        end
        names = strsplit(m{mid},'_');
        YDat(1,(mid-1)*takenrt+1:mid*takenrt) = repmat(find(strcmp(BGOptNames,names{1})),[takenrt,1]);
        YDat(2,(mid-1)*takenrt+1:mid*takenrt) = repmat(find(strcmp(SideOpt,names(2))),[takenrt,1]);
    end
    
    XDat(:,sum(isnan(XDat),1)>0) = []; %REmove any pixels that have missing data or have too low variance    
    Results = GP_MTL_CLASS(XDat,YDat,0,nfolds,randomtidx,2,0);
    
    
    %Save Results
    save(fullfile(ResultPath, [mouse 'time' num2str(TW{TWidx}(1)) '-' num2str(TW{TWidx}(2)), Model.Rnames{regio2take(roidx)}]),'Results')
    
    
end
end
