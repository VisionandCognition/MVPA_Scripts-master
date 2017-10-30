function SVMWideField_StimDecoder(Mouseidx,TWidx,ridx)
% Decoding --> Use MVPA (SVM) to decode orientation of the stimulus for one
% session
if ~exist('Mouseidx','var')
    Mouseidx =4;
    TWidx = 2;
end
nboot=100
wholebrainana = 1
spotlightana = 1
nrSVM = 12;
if ischar(Mouseidx)
    Mouseidx = str2num(Mouseidx);
end
if ischar(TWidx)
    TWidx = str2num(TWidx);
end
xpix =800;
ypix = 800;
%Station
clear notused
%% User Defined Input
miceopt = {'Alladin','Chief','Esmeralda','Frey'} %options for mice
Stim2Check = 'DelayedOriTuningSound'%Name of the stimulus as written in the LOG-file
%Timelimit: Don't need data from time after this.
timelimit1 = 2500; %ms
ReactionOpt = {'Error','Hit','Miss'};
if strcmp(Stim2Check,'DelayedOriTuningSound')
    basel = [-300,-50];
    TW = {[-300,-50],[100,500],[750,1350],[1350,1950]};
    TWNames = {'Baseline','Visual','Delay','Response'};
else
    basel = [-250 -50];
    fgtw = [120 250]; %FOR FG
    VisInit = [50 120];
    bigtw = [200 450];
    TW = {basel,VisInit,fgtw,bigtw}
    TWNames = {'Baseline','VisualInit','fgtw','bigtw'};
end
cd
if ispc
    DataPath = '\\vcnin\mouse_working_memory\Data4Class\TMPData'; % Set path
    ResultPath = 'I:\SARA\TMPResults'; % Set path
    ScriptsPath = 'I:\SARA\MVPA_Scripts'; %Set path
    TWidx = 3;
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

tempstation = fullfile(ResultPath,'TMPMatlab');
if ~exist(tempstation,'dir')
    mkdir(tempstation)
end
%Load info
load(fullfile(DataPath,'sessionstruct.mat'))
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
    %Load Alan Brain model
    BrainModel{midx} = load(fullfile(DataPath,'brainareamodel.mat'))
    for ridx = ridx
        PERF = cell(1,nboot+1);
        SVMMAP = PERF;
        notzeroanymoreright = 0;
        notzeroanymoreleft = 0;

        %Create a file on the disk
        delete(fullfile(tempstation,[mouse  '_TW' num2str(TWidx) ReactionOpt{ridx} '_tmpfile.mat']))
        TMPMAT = matfile(fullfile(tempstation,[mouse '_TW' num2str(TWidx) ReactionOpt{ridx} '_tmpfile.mat'])); %Make a workable matfile
       
        sessioncount = 0;
        clear RawData
        rightdat = [];
        leftdat = [];
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
                    load(fullfile('\\vcnin\mouse_working_memory\Imaging\',mouse,[mouse date],[mouse num2str(expnr)],[mouse num2str(expnr) '.mat']));
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
                        CheckReactions([this.folder this.expname '\' this.mouse this.expnum '.mat'])
                        tmp = load([folder expname '\' mouse expnum '.mat']);
                        if isfield(tmp,'tosave')
                            tmp = tmp.tosave;
                        end
                    end
                    LOG.Reaction = LOG.correctReaction; %Change the reactions into checked reactions
                end
                
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
                ConditionNames = ConditionNames(cvec);
                
                idx = find(~cellfun(@isempty,strfind(ConditionNames,'Too Early')));
                ConditionNames(idx) = cellfun(@(X) strrep(X,X(strfind(X,'Too Early'):9),'TooEarly'),ConditionNames(idx),'UniformOutput',0);
                
                %Average over orientations
                conditionparts = cellfun(@(X) strsplit(X,' '),ConditionNames,'UniformOutput',0);
                
                %Find all reactions
                reaction = cellfun(@(X) X{1},conditionparts,'UniformOutput',0); %Reaction
                orientation = cellfun(@(X) X{2},conditionparts,'UniformOutput',0); %orientations
                OriOpt = unique(orientation);
                side = cellfun(@(X) X{4},conditionparts,'UniformOutput',0); %SIdes
                SideOpt = unique(side);
                
                %% Load data movement matrix
                load(fullfile(DataPath,[mouse date],[mouse num2str(expnr)],'ThrowAwayIdx.mat'))
                if strcmp(Stim2Check,'DelayedOriTuningSound')
                    fullfgtr = find(LOG.currentdelay==1500 & LOG.Gavepassive(LOG.RealTask==1) == 0 & LOG.Ezbox == 0& LOG.TotalStimDur == 500);
                else
                    if strcmp(trialtypes{id},'FG')
                        fullfgtr = find(LOG.BGContrast==1 & LOG.Gavepassive==0&LOG.Ezbox==0 & LOG.OOP ==0);
                    elseif strcmp(trialtypes{id},'GREY')
                         fullfgtr = find(LOG.BGContrast==0 & LOG.Gavepassive==0&LOG.Ezbox==0);
                    elseif strcmp(trialtypes{id},'OOP')
                        fullfgtr = find(LOG.BGContrast==1 & LOG.Gavepassive==0&LOG.Ezbox==0 & LOG.OOP ==1);

                    end
                end
                
                %%  Load 'drift correction'
                load(fullfile(DataPath,[mouse date],[mouse num2str(expnr)],'BASELINEMAT.mat'))
                
                %% Load data for timeline etc.
                rawdatfiles = dir(fullfile(DataPath,[mouse date],[mouse num2str(expnr)],[mouse num2str(expnr) '_RawData*']));
                
                load(fullfile(DataPath,[mouse date],[mouse num2str(expnr)],rawdatfiles(strcmp({rawdatfiles(:).name},[mouse num2str(expnr) '_RawData_C' num2str(length(cvec)) '.mat'])).name));
                clear conddata
                twidx = find(timeline>=TW{TWidx}(1)&timeline<=TW{TWidx}(2));
                baseidx = find(timeline>=basel(1) & timeline<=basel(2));
                
                %% load in rawdata
                
                %% RIGHT
                rightidx = find(~cellfun(@isempty,(cellfun(@(X) strfind(X,'right'),ConditionNames,'UniformOutput',0)))& ~cellfun(@isempty,(cellfun(@(X) strfind(X,ReactionOpt{ridx}),ConditionNames,'UniformOutput',0))));
                if strcmp(Stim2Check,'DelayedOriTuningSound')
                    rightidx(ismember(rightidx,find(~cellfun(@isempty,(cellfun(@(X) strfind(X,'500'),ConditionNames,'UniformOutput',0)))))) = [];
                end
                for i = 1:length(rightidx)
                    tmpload = load(fullfile(DataPath,[mouse date],[mouse num2str(expnr)],rawdatfiles(strcmp({rawdatfiles(:).name},[mouse num2str(expnr) '_RawData_C' num2str(rightidx(i)) '.mat'])).name));
                    tmpload.conddata = single(tmpload.conddata); %uint8 == 0 means single =nan;
                    tmpload.conddata(tmpload.conddata==0)=nan;
                    try
%                         rmtmp = ~removeidx(1:length(tmpload.ctrials{rightidx(i)}),rightidx(i))';
%                         rm2tmp = ismember(tmpload.ctrials{rightidx(i)},fullfgtr);
                        rm3tmp = ismember(tmpload.ctrials{rightidx(i)},fullfgtr);
%                         do this one later
%                         rm3tmp = (rmtmp==1 & rm2tmp==1); %Do this one
%                         later
                        trialidx = tmpload.ctrials{rightidx(i)};
                        trialidx = trialidx(rm3tmp);
                        rightdattmp =  zeros(800,800,length(twidx),sum(rm3tmp),'single');
                        for j = 1:100:xpix
                            tmpnw = single(tmpload.conddata(j:j+99,:,:,rm3tmp))./permute(repmat(BASELINEMAT(j:j+99,:,trialidx),[1,1,1,size(tmpload.conddata,3)]),[1,2,4,3]);
                            rightdattmp(j:j+99,:,:,:) = (tmpnw(:,:,twidx,:)-repmat(nanmean(tmpnw(:,:,baseidx,:),3),[1,1,length(twidx),1]))./repmat(nanmean(tmpnw(:,:,baseidx,:),3),[1,1,length(twidx),1]);
                        end
                        
                        %Remove trials with black windows
                        rightdattmp = rightdattmp(:,:,:,~(sum(isnan(reshape(nanmean(rightdattmp,3),[xpix*ypix,size(rightdattmp,4)])),1)>0.7*(xpix*ypix)));
                        %Remove trials with inf average values (can happen when
                        %screen was black?)
                        tmpf = nanmean(reshape(rightdattmp,[size(rightdattmp,1)*size(rightdattmp,2)*size(rightdattmp,3),size(rightdattmp,4)]));
                        rightdattmp(:,:,:,tmpf>0.5) = [];
                        trialidx(tmpf>0.5) = [];
                        
                        if isempty(trialidx)
                            warning('No trials for right in this session left')
                            nrtotalRight{sessioncount} = 0;
                            continue
                        else
                           notzeroanymoreright = 1; 
                        end
                        if (sessioncount==1 && i == 1) || ~notzeroanymoreright
                            nrtotalRight{sessioncount} = length(trialidx);
                        elseif i ~= 1
                            nrtotalRight{sessioncount} = nrtotalRight{sessioncount}+length(trialidx);
                        else
                            nrtotalRight{sessioncount} = nrtotalRight{sessioncount-1}+length(trialidx);
                        end                                                
                        
                        m = whos(TMPMAT);
                        m = {m(:).name};
                        
                        if ~ismember('rightdat',m)
                            if length(size(rightdattmp))<4
                                rightdattmp = cat(4,rightdattmp,nan(size(rightdattmp)));
                            end
                           TMPMAT.rightdat = squeeze(nanmean(rightdattmp,3));
                        elseif i ~= 1
                            TMPMAT.rightdat(:,:,nrtotalRight{sessioncount} - sum(trialidx)+1:nrtotalRight{sessioncount}) = squeeze(nanmean(rightdattmp,3));
                        else
                            TMPMAT.rightdat(:,:,nrtotalRight{sessioncount-1}+1:nrtotalRight{sessioncount}) = squeeze(nanmean(rightdattmp,3));                            
                        end
                        
                    catch ME
                      
                        disp(ME)
                        if strcmp(ME.identifier,'MATLAB:nomem')
                           continue 
                        else
                        keyboard
                        end
                    end
                    clear tmpload
                    clear rightdattmp
                end
               
                
                %% LEFT
                leftidx = find(~cellfun(@isempty,(cellfun(@(X) strfind(X,'left'),ConditionNames,'UniformOutput',0))) & ~cellfun(@isempty,(cellfun(@(X) strfind(X,ReactionOpt{ridx}),ConditionNames,'UniformOutput',0))));
                if strcmp(Stim2Check,'DelayedOriTuningSound')
                    leftidx(ismember(leftidx,find(~cellfun(@isempty,(cellfun(@(X) strfind(X,'500'),ConditionNames,'UniformOutput',0)))))) = [];
                end
                for i = 1:length(leftidx)
                    tmpload = load(fullfile(DataPath,[mouse date],[mouse num2str(expnr)],rawdatfiles(strcmp({rawdatfiles(:).name},[mouse num2str(expnr) '_RawData_C' num2str(leftidx(i)) '.mat'])).name));
                    tmpload.conddata = single(tmpload.conddata);
                    tmpload.conddata(tmpload.conddata==0)=nan;
                    try
%                         rmtmp = ~removeidx(1:length(tmpload.ctrials{leftidx(i)}),leftidx(i))';
%                         rm2tmp = ismember(tmpload.ctrials{leftidx(i)},fullfgtr);;
%                         rm3tmp = (rmtmp==1 & rm2tmp==1); %Do this one
% % % %                         later
                        rm3tmp = ismember(tmpload.ctrials{leftidx(i)},fullfgtr); %Do this one later
                        
                        trialidx = tmpload.ctrials{leftidx(i)};
                        trialidx = trialidx(rm3tmp);
                        leftdattmp = nan(800,800,length(twidx),sum(rm3tmp),'single');
                        for j = 1:100:xpix
                            tmpnw =  single(tmpload.conddata(j:j+99,:,:,rm3tmp))./permute(repmat(BASELINEMAT(j:j+99,:,trialidx),[1,1,1,size(tmpload.conddata,3)]),[1,2,4,3]);
                            leftdattmp(j:j+99,:,:,:) = (tmpnw(:,:,twidx,:)-repmat(nanmean(tmpnw(:,:,baseidx,:),3),[1,1,length(twidx),1]))./repmat(nanmean(tmpnw(:,:,baseidx,:),3),[1,1,length(twidx),1]);
                        end
                        %Remove trials with black windows
                        leftdattmp = leftdattmp(:,:,:,~(sum(isnan(reshape(nanmean(leftdattmp,3),[xpix*ypix,size(leftdattmp,4)])),1)>0.7*(xpix*ypix)));
                        
                        %Remove trials with inf average values (can happen when
                        %screen was black?)
                        tmpf = nanmean(reshape(leftdattmp,[size(leftdattmp,1)*size(leftdattmp,2)*size(leftdattmp,3),size(leftdattmp,4)]));
                        leftdattmp(:,:,:,tmpf>0.5) = [];
                        
                        trialidx(tmpf>0.5) = [];
                          if isempty(trialidx)
                            warning('No trials for left in this session left')
                             nrtotalLeft{sessioncount} = 0;
                            continue
                        else
                           notzeroanymoreleft = 1; 
                        end
                        if (sessioncount==1 && i == 1) || ~notzeroanymoreleft                        
                           nrtotalLeft{sessioncount} = length(trialidx);
                        elseif i ~= 1
                            nrtotalLeft{sessioncount} = nrtotalLeft{sessioncount}+length(trialidx);
                        else
                            nrtotalLeft{sessioncount} = nrtotalLeft{sessioncount-1}+length(trialidx);
                        end                                                
                        
                        m = whos(TMPMAT);
                        m = {m(:).name};
                        
                        if ~ismember('leftdat',m)
                            if length(size(leftdattmp))<4
                                leftdattmp = cat(4,leftdattmp,nan(size(leftdattmp)));
                            end
                           TMPMAT.leftdat = squeeze(nanmean(leftdattmp,3));
                        elseif i ~= 1
                            TMPMAT.leftdat(:,:,nrtotalLeft{sessioncount} - sum(trialidx)+1:nrtotalLeft{sessioncount}) = squeeze(nanmean(leftdattmp,3));
                        else
                            TMPMAT.leftdat(:,:,nrtotalLeft{sessioncount-1}+1:nrtotalLeft{sessioncount}) = squeeze(nanmean(leftdattmp,3));                            
                        end                        
                    catch ME
                        if strcmp(ME.identifier,'MATLAB:nomem')
                            continue
                        else
%                             keyboard
                        end
                    end
                    clear tmpload
                    clear leftdattmp
                end
           
            end
        end
        
        %% SVM
        rightdat = TMPMAT.rightdat;
        leftdat = TMPMAT.leftdat;
        
        delete(fullfile(tempstation,[mouse  '_TW' num2str(TWidx) ReactionOpt{ridx} '_tmpfile.mat']))
       
          
        trainp = 0.9;
        nleft = size(leftdat,3);
        nright =  size(rightdat,3);
        
        if nleft>nright
            tn = floor(trainp.*nright); %Number of trial in training set
            tt = nright-tn; %Number of trials in test set
        else
            tn = floor(trainp.*nleft); %Number of trial in training set
            tt = nleft-tn; %Number of trials in test set
        end
        
        %Remove trials because its getting'too big to handle
        if tn > 150
            disp('More than 150 trials')
            tn = 150;
            tt = 0.1*tn;
            randl = randperm(nleft,tn+tt);
            randr = randperm(nright,tn+tt);
            
            rightdat = rightdat(:,:,randr);
            nright = size(rightdat,3);
            
            leftdat = leftdat(:,:,randl);
            nleft = size(leftdat,3);
        end
        
        %Remove areas
       %Remove areas
        throwawayareas = find(cellfun(@isempty,BrainModel{midx}.Model.Rnames));
        throwawayareas = [throwawayareas; find(cellfun(@(X) ismember(X,{'OlfactoryBulb','fibrtracts','InfCol','SupColSens'}),BrainModel{midx}.Model.Rnames))];
        keepareas = 1:length(BrainModel{midx}.Model.Rnames);
        keepareas(throwawayareas)=[];
        removepix = true(xpix,ypix);
        for areaid = 1:length(keepareas)
            bounds = BrainModel{midx}.Model.Boundaries{keepareas(areaid)};
            for boundid = 1:length(bounds)
                removepix(poly2mask(bounds{boundid}(:,1),bounds{boundid}(:,2),xpix,ypix)) = 0;
            end
        end
        
        removepix = smooth2a(removepix,5);
        removepix(removepix<0.9)=0;
        removepix = ~imfill(~removepix,'holes');

        removepixvec = reshape(removepix,[xpix*ypix,1]);
        if wholebrainana
            thistimer = tic;
            tmpright = reshape(rightdat,[xpix*ypix,nright]);
            tmpleft = reshape(leftdat,[xpix*ypix,nleft]);
            
            tmpleft = tmpleft';
            tmpright = tmpright';
            
            removenanpix = find(squeeze(sum(isnan(tmpleft),1)>0) | squeeze(sum(isnan(tmpright),1)>0) | removepixvec' == 1);
            tmpleft(:,removenanpix) = [];
            tmpright(:,removenanpix) = [];
            
            if isempty(tmpleft) || isempty(tmpright)
                continue
            end
            
            
            %% Whole-brain approach
            %Run svm
            clear SVMModel
            beta = zeros(xpix*ypix-length(removenanpix),nrSVM);
            %24 SVMs (runs fast on 12 cores)
            svmperf = nan(1,nrSVM);
            parfor s = 1:nrSVM
                %Left - select the trials that go into the errors
                leftperm = randperm(nleft);
                leftpicktr = leftperm(1:tn);
                leftpickte = leftperm(tn+1:tn+tt); %Check this line
                
                lefttrain = tmpleft(leftpicktr,:);%training data
                lefttest =  tmpleft(leftpickte,:);%test data
                %Right trials
                rightperm = randperm(nright);
                rightpicktr = rightperm(1:tn);
                rightpickte = rightperm(tn+1:tn+tt); %Check this line
                righttrain =  tmpright(rightpicktr,:); %training
                righttest = tmpright(rightpickte,:); %test
                
                %Train SVM with linear kernel, standradise (z-score) predictors
                trainset = [lefttrain;righttrain];
                out = [ones(tn,1);ones(tn,1).*-1];
                try
                    
                    SVMModel = fitcsvm(trainset,out,'Standardize','on','KernelFunction','linear','KernelScale','auto','BoxConstraint',1);
                    
                    %Now test the model with the left over trials
                    testset = [lefttest;righttest];
                    class = [ones(tt,1);ones(tt,1).*-1];
                    label = predict(SVMModel,testset);
                    svmperf(s) = mean(class==label); %Performance of SVM
                    beta(:,s) = SVMModel.Beta; %Feature weights
                catch ME
                    disp(ME)
                end
            end
            
            %The absoulute feature weights
            betaform = true(xpix*ypix,1);
            betaform(removenanpix') = 0;
            newbeta = nan(xpix*ypix,1);
            newbeta(betaform) = smooth2a(mean((beta),2),2);% abs
            
            PERF{1}.svmperf = svmperf;
            PERF{1}.beta = reshape(newbeta,xpix,ypix);
            disp(['Original Estimate whole-brain approach for ' num2str(TW{TWidx}(1)) '-' num2str(TW{TWidx}(2)) ' took ' num2str(toc(thistimer)./60) ' minutes'])
            
            clear SVMModel
            clear beta
            %% Now we bootstrap (resample from the same data with replacement to form the 0-distribution)
            rng default
            parfor bi = 1:nboot
                beta = zeros(xpix*ypix-length(removenanpix),nrSVM);
                %24 SVMs (runs fast on 12 cores)
                svmperf = nan(1,nrSVM);
                for s = 1:nrSVM
                    %Left - select the trials that go into the errors
                    leftperm = randi(nleft,1,nleft);
                    leftpicktr = leftperm(1:tn);
                    leftpickte = leftperm(tn+1:tn+tt); %Check this line
                    
                    lefttrain = tmpleft(leftpicktr,:);%training data
                    lefttest =  tmpleft(leftpickte,:);%test data
                    
                    %Right trials
                    rightperm = randi(nright,1,nright);
                    rightpicktr = rightperm(1:tn);
                    rightpickte = rightperm(tn+1:tn+tt); %Check this line
                    righttrain =  tmpright(rightpicktr,:); %training
                    righttest = tmpright(rightpickte,:); %test
                    
                    %Train SVM with linear kernel, standradise (z-score) predictors
                    trainset = [lefttrain;righttrain];
                    out = [ones(tn,1);ones(tn,1).*-1];
                    try
                        
                        SVMModel = fitcsvm(trainset,out,'Standardize','on','KernelFunction','linear','KernelScale','auto','BoxConstraint',1);
                        
                        %Now test the model with the left over trials
                        testset = [lefttest;righttest];
                        class = [ones(tt,1);ones(tt,1).*-1];
                        label = predict(SVMModel,testset);
                        svmperf(s) = mean(class==label); %Performance of SVM
                        beta(:,s) = SVMModel.Beta; %Feature weights
                    catch ME
                        disp(ME)
                    end
                end
                %The absoulute feature weights
                newbeta = nan(xpix*ypix,1);
                newbeta(betaform) = smooth2a(mean((beta),2),2);% abs
                
                PERF{1+bi}.svmperf = svmperf;
                PERF{1+bi}.beta =  reshape(newbeta,xpix,ypix);
               
            end
          
            clear SVMModel
             %Save Results
             save(fullfile(ResultPath, [mouse '_' ReactionOpt{ridx} 'time' num2str(TW{TWidx}(1)) '-' num2str(TW{TWidx}(2))]),'PERF','SVMMAP','miceopt','ReactionOpt','TW','tt','tn','-v7.3')

        end
         if spotlightana
            %% Use a spotlighht based approach
            squaresize = 25;
            colsize = ypix;
            rowsize = xpix;
            colstart = 1:squaresize/5:colsize-squaresize;
            rowstart = 1:squaresize/5:rowsize-squaresize;
            [coli,rowi] = meshgrid(colstart,rowstart);
            svmmap = zeros(length(rowstart),length(colstart));
            thistimer = tic;
            for i = 1:length(colstart)
                for j = 1:length(rowstart)
                    tmpright = reshape(rightdat(rowstart(j):rowstart(j)+squaresize-1,colstart(i):colstart(i)+squaresize-1,:),[squaresize*squaresize,nright]);
                    tmpright = tmpright';
                    tmpleft = reshape(leftdat(rowstart(j):rowstart(j)+squaresize-1,colstart(i):colstart(i)+squaresize-1,:,:),[squaresize*squaresize,nleft]);
                    tmpleft = tmpleft';
                    removenanpix = find(squeeze(sum(isnan(tmpleft),1)>0) | squeeze(sum(isnan(tmpright),1)>0) | reshape(removepix(rowstart(j):rowstart(j)+squaresize-1,colstart(i):colstart(i)+squaresize-1),[squaresize*squaresize,1])' == 1);
                    tmpleft(:,removenanpix) = [];
                    tmpright(:,removenanpix) = [];
                    if isempty(tmpleft) || isempty(tmpright)
                        continue
                    end                    
                    svmperf = nan(1,nrSVM);
                    parfor s = 1:nrSVM
                        %Left - select the trials that go into the errors
                        leftperm = randperm(nleft);
                        leftpicktr = leftperm(1:tn);
                        leftpickte = leftperm(tn+1:tn+tt); %Check this line
                        
                        lefttrain = tmpleft(leftpicktr,:);%training data
                        lefttest =  tmpleft(leftpickte,:);%test data
                        %Right trials
                        rightperm = randperm(nright);
                        rightpicktr = rightperm(1:tn);
                        rightpickte = rightperm(tn+1:tn+tt); %Check this line
                        righttrain =  tmpright(rightpicktr,:); %training
                        righttest = tmpright(rightpickte,:); %test
                        
                        %Train SVM with linear kernel
                        trainset = [lefttrain;righttrain];
                        out = [ones(tn,1);ones(tn,1).*-1];
                        SVMModel = fitcsvm(trainset,out,'Standardize','on','KernelFunction','linear','KernelScale','auto','BoxConstraint',1);
                        
                        %Now test the model with the left over trials
                        testset = [lefttest;righttest];
                        class = [ones(tt,1);ones(tt,1).*-1];
                        label = predict(SVMModel,testset);
                        svmperf(s) = mean(class==label);
                    end
                    
                    svmmap(j,i) = mean(svmperf);
                end
            end
            try
                SVMMAP{1} = svmmap;
            catch ME
                disp(ME)
%                 keyboard
            end
            disp(['Spotlight approach for ' num2str(TW{TWidx}(1)) '-' num2str(TW{TWidx}(2))  ' took ' num2str(toc(thistimer)./60) ' minutes'])
            
            
            clear SVMModel
            clear beta
            clear svmmap
            %% Now we bootstrap (resample from the same data with replacement to form the 0-distribution)
            rng default
            try
                parfor bi = 1:nboot                    
                    squaresize = 25;
                    colsize = ypix;
                    rowsize = xpix;
                    colstart = 1:squaresize/5:colsize-squaresize;
                    rowstart = 1:squaresize/5:rowsize-squaresize;
                    [coli,rowi] = meshgrid(colstart,rowstart);
                    svmmap = zeros(length(rowstart),length(colstart));
                    thistimer = tic;
                    for i = 1:length(colstart)
                        for j = 1:length(rowstart)
                            tmpright = reshape(rightdat(rowstart(j):rowstart(j)+squaresize-1,colstart(i):colstart(i)+squaresize-1,:),[squaresize*squaresize,nright]);
                            tmpright = tmpright';
                            tmpleft = reshape(leftdat(rowstart(j):rowstart(j)+squaresize-1,colstart(i):colstart(i)+squaresize-1,:,:),[squaresize*squaresize,nleft]);
                            tmpleft = tmpleft';
                            removenanpix = find(squeeze(sum(isnan(tmpleft),1)>0) | squeeze(sum(isnan(tmpright),1)>0) | reshape(removepix(rowstart(j):rowstart(j)+squaresize-1,colstart(i):colstart(i)+squaresize-1),[squaresize*squaresize,1])' == 1);
                            tmpleft(:,removenanpix) = [];
                            tmpright(:,removenanpix) = [];
                            if isempty(tmpleft) || isempty(tmpright)
                                continue
                            end
                            
                            svmperf = nan(1,nrSVM);
                            for s = 1:nrSVM
                                %Left - select the trials that go into the errors
                                leftperm = randi(nleft,1,nleft);
                                leftpicktr = leftperm(1:tn);
                                leftpickte = leftperm(tn+1:tn+tt); %Check this line
                                
                                lefttrain = tmpleft(leftpicktr,:);%training data
                                lefttest =  tmpleft(leftpickte,:);%test data
                                %Right trials
                                rightperm = randi(nright,1,nright);
                                rightpicktr = rightperm(1:tn);
                                rightpickte = rightperm(tn+1:tn+tt); %Check this line
                                righttrain =  tmpright(rightpicktr,:); %training
                                righttest = tmpright(rightpickte,:); %test
                                
                                %Train SVM with linear kernel
                                trainset = [lefttrain;righttrain];
                                out = [ones(tn,1);ones(tn,1).*-1];
                                SVMModel = fitcsvm(trainset,out,'Standardize','on','KernelFunction','linear','KernelScale','auto','BoxConstraint',1);
                                
                                %Now test the model with the left over trials
                                testset = [lefttest;righttest];
                                class = [ones(tt,1);ones(tt,1).*-1];
                                label = predict(SVMModel,testset);
                                svmperf(s) = mean(class==label);
                            end
                            
                            svmmap(j,i) = mean(svmperf);
                        end
                    end
                    
                    SVMMAP{1+bi} = svmmap;
                end
            catch ME
                disp(ME)
%                 keyboard
            end
            
         end
        %Save Results
        save(fullfile(ResultPath, [mouse '_' ReactionOpt{ridx} 'time' num2str(TW{TWidx}(1)) '-' num2str(TW{TWidx}(2))]),'PERF','SVMMAP','miceopt','ReactionOpt','TW','tt','tn','-v7.3')
    
    end 
end
end
