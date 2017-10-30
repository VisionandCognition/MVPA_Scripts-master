function FG_MVPA_SVM_RFE(Mouseidx,TWidx,roidx,roimethod)
% Decoding --> Use MVPA (SVM) to decode orientation of the stimulus for one
% session
if ~exist('Mouseidx','var')
    Mouseidx =1;
    TWidx = 2;
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
miceopt = {'Jules','Marsellus','Vincent','Zed'} %options for mice
Stim2Check = 'FGTask'%Name of the stimulus as written in the LOG-file
%Timelimit: Don't need data from time after this.
timelimit1 = 2500; %ms
BGOpt = [0,1];
taskOptNames = {'FG','Side'};
BGOptNames = {'Grey','Contrast'};
nrtask = 2;
if strcmp(Stim2Check,'FGTask')
    basetw = {[-300 -100]};
    TW = {[-100,50],[100,250],[250 500],[500 750]};
    TWNames = {'Baseline','Visual'};
end
cd
if ispc
    DataPath = '\\vcnin\mouse_working_memory\Data4Class\TMPData'; % Set path
    ResultPath = 'I:\SARA\TMPResults'; % Set path
    ScriptsPath = 'I:\SARA\MVPA_Scripts'; %Set path
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

addpath(fullfile(ScriptsPath,'\libsvm-3.20\matlab\'))

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
        
        %Cut out brain region
        mask = zeros(800,800);
        for entirebrainid = 11:length(Model.Boundaries)-2
            Borders = Model.Boundaries{entirebrainid};
            for roi2dx = 1:length(Borders)
                masktmp = poly2mask(Borders{roi2dx}(:,1),Borders{roi2dx}(:,2),800,800);
                %Shrink to not have border effects
                masktmp = bwmorph(masktmp,'shrink',1);
                mask(masktmp==1)=1;
            end
        end
    end
    
    if roidx > length(regio2take)
        disp(['too large number for the number of rois that exist'])
        return
    end
    if roidx == length(regio2take)
        disp(['Check whether all rois are included..'])
    end
    
    
    delete(fullfile(tempstation,[mouse  '_TW' num2str(TWidx) '_tmpfile.mat']))
    TMPMAT = matfile(fullfile(tempstation,[mouse '_TW' num2str(TWidx) '_tmpfile.mat'])); %Make a workable matfile
    
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
                            
                            base = squeeze(single(nanmean(conddata(:,:,(timeline>=basetw{1}(1) & timeline<=basetw{1}(2)),keepidx),3)));
                            tmp = single(conddata(:,:,(timeline>=TW{TWidx}(1) & timeline <= TW{TWidx}(2)),keepidx))./permute(repmat(base,[1,1,1,sum((timeline>=TW{TWidx}(1) & timeline <= TW{TWidx}(2)))]),[1,2,4,3]); %Shiftdim!
                            tmp(tmp==0)=nan;                                                        
                            tmp = squeeze(nanmean(tmp,3));
                            tmp(~repmat(mask,[1,1,size(tmp,3)])) = nan;    

                            tmpname = [BGOptNames{ridx} '_' SideOpt{stidx}];
                            if ~ismember(tmpname,m)
                                if length(size(tmp))<3
                                    tmp = cat(3,tmp,nan(size(tmp)));
                                end
                                eval(['TMPMAT.' tmpname ' = tmp;']);
                                clear tmp
                            elseif cidxcount ~= 1
                                    eval(['TMPMAT.' tmpname '(:,:,nrtotal{sessioncount,stidx,ridx}-sum(keepidx)+1:nrtotal{sessioncount,stidx,ridx}) = tmp;'])
                            else
                                    eval(['TMPMAT.' tmpname '(:,:,nrtotal{sessioncount-1,stidx,ridx}+1:nrtotal{sessioncount,stidx,ridx}) = tmp;'])
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
            
            if sessioncount >= 1
                break
            end
        end
        if sessioncount >= 1
            break
        end
    end
    clear conddata
    clear SErrorval
    clear RawTmp
    
    nrtrials = nrtotal(sessioncount,:,1);
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
    shamtrialvec1 = randperm(takenrt);
    shamtrialvec2 = randperm(takenrt);
    %% Start Decoding
    % Add paths necessary
    %Call TMPMAT
    TMPMAT = matfile(fullfile(tempstation,[mouse '_TW' num2str(TWidx) '_tmpfile.mat'])); %Make a workable matfile
    
    % Write data into GP understandable 'structs'
    m = whos(TMPMAT);
    m = {m(:).name};
    nrtask = 1;

    xpix = size(eval(['TMPMAT.' m{1}]),1);
    ypix = size(eval(['TMPMAT.' m{1}]),2);
    try
        ReactionNumbers = 1:length(unique(BGOpt));
        SideNumbers = 1:length(unique(SideOpt));
        p=4;
        XDat = zeros((length(m)-2)*takenrt,(xpix*ypix)/(p^2),'single');
        YDat = zeros(nrtask,(length(m)-2)*takenrt,'single');
        for mid = 1:length(m)
            if ~isempty(strfind(m{mid},'Contrast'))
                continue
            end
            tmp =  eval(['TMPMAT.' m{mid}]);
            randomtake = shamtrialvec(find(shamtrialvec<=size(tmp,3)));
            tmp = tmp(:,:,randomtake(1:takenrt));
            
            % Downsample
            for trialidx = 1:takenrt
                M = tmp(:,:,trialidx);
                [xpix,ypix]=size( M); %M is the original matrix
                
                M= sum(reshape( M,p,[]) ,1 );
                M=reshape(M,xpix/p,[]).'; %Note transpose
                M=sum(reshape(M,p,[]) ,1);
                tmp2 = reshape(reshape(M,ypix/p,[]),[xpix/4*ypix/4,1]);
                XDat((mid-3)*takenrt+trialidx,:) = tmp2;
            end
            names = strsplit(m{mid},'_');
%             YDat(1,(mid-1)*takenrt+1:mid*takenrt) = repmat(find(strcmp(BGOptNames,names{1})),[takenrt,1]);
            YDat(1,(mid-3)*takenrt+1:(mid-2)*takenrt) = repmat(find(strcmp(SideOpt,names(2))),[takenrt,1]);
        end
        throwpix = find(sum(isnan(XDat),1)>0);
        XDat(:,sum(isnan(XDat),1)>0) = []; %REmove any pixels that have missing data or have too low variance
        
        %Z-score
        stop = 0;
        count = 1;
        while count<10 && ~stop
            ZSc = (XDat - repmat(nanmean(XDat,1),[size(XDat,1),1]))./repmat(nanstd(XDat,[],1),[size(XDat,1),1]);
            if any(abs(ZSc(:))>2)
                ZSc(ZSc>2)= 2;
                ZSc(ZSc<-2)=-2;
            else
                stop=1;
            end
            count = count+1;
            XDat = ZSc;

        end
        
        %Remove trials that are constantly at max
        rmtrial = find((sum(XDat>=2,2)>=0.9*size(XDat,1))|sum(XDat<=-2,2)>=0.9*size(XDat,1));
        XDat(rmtrial,:) = [];
        YDat(rmtrial) = [];
        
         %  Take 50% most variant pixels
         varorder = sort(var(XDat,[],1),'ascend');
        rmpixel = find(var(XDat,[],1) > varorder(floor(length(varorder)./(2*3))));
        XDat(:,rmpixel) = [];
        randomtidx(find(randomtidx > size(XDat,1))) = [];
        chunksz = floor(length(randomtidx)./nfolds); %Even chunks of data as nested-cross

    catch ME
        disp(ME)
        keyboard
    end
    actvspred = nan(nrtask,length(randomtidx),2);
    %% Actual SVM + Cross-val
    TOOKC = nan(1,nfolds);
    for nc = 1:nfolds
        testcidx = randomtidx((nc-1)*chunksz+1:nc*chunksz); %fold-testdata
        trainingcidx = randomtidx(~ismember(randomtidx,testcidx)); %fold-trainingdata
        
        training = XDat(trainingcidx,:);
        maxdat = nanmax(training);
        mindat = nanmin(training);
        % Normalize between 0 and 1
        for i = 1:100:size(training,1)
            try
                training(i:i+99,:) = (training(i:i+99,:)-repmat(mindat,[100,1]))./repmat((maxdat-mindat),[100,1]);
            catch
                training(i:end,:) = (training(i:end,:) - repmat(mindat,[size(training,1)-i+1,1]))./repmat((maxdat-mindat),[size(training,1)-i+1,1]);
            end
        end
        trainlab = YDat(:,trainingcidx); %label
        
        figure; subplot(2,2,1); imagesc(training(trainlab==1,:))
        title('left')
        subplot(2,2,2); imagesc(training(trainlab==2,:))
        title('right')
        ylabel(['Mouse ' mouse '_TW' num2str(TWidx) 'Fold ' num2str(nc)])
        %Labels for this fold (select part of trainingset as valid.set
        test = XDat(testcidx,:);
        % Normalize between 0 and 1
        test = (test-repmat(mindat,[size(test,1),1]))./repmat((maxdat-mindat),[size(test,1),1]);
        testlab = YDat(:,testcidx);
        
        pixelmap = nan(xpix/p,ypix/p);
        keeppix = 1:prod(size(pixelmap));
        keeppix(throwpix) = [];
        keeppix(rmpixel) = [];
        pixelmap(keeppix)= nanmean(training(trainlab==1,:),1);
        subplot(2,2,3); imagesc(pixelmap')
        
        pixelmap = nan(xpix/p,ypix/p);
        keeppix = 1:prod(size(pixelmap));
        keeppix(throwpix) = [];
        keeppix(rmpixel) = [];

        pixelmap(keeppix)= nanmean(training(trainlab==2,:),1);
        subplot(2,2,4); imagesc(pixelmap') 
        
        try
            %% Optimize parameter c on sub-set
            accuracy = nan(nrtask,length(n));
            mse = accuracy;
            figure;
            nestedtestidx = randperm(size(trainlab,2),floor(size(trainlab,2)./10));
            nestedtrainidx = 1:size(trainlab,2);
            nestedtrainidx(nestedtestidx) = [];
            for t = 1:nrtask
                for i = 1:numel(n) %
                    c = 2^n(i);
                    model = svmtrain(double(trainlab(t,nestedtrainidx))',double(training(nestedtrainidx,:)),['-q -c ' num2str(c)]);
                    [lbl, acc, dec] = svmpredict(double(trainlab(t,nestedtestidx))',double(training(nestedtestidx,:)),model,['-q']);
                    accuracy(t,i) = acc(1);
                    mse(t,i) = acc(2);
                end
                [acc_maxid] = find(accuracy(t,:)==nanmax(accuracy(t,:)));
                [MSE_minid] = find(mse(t,:)==nanmin(mse(t,:)));
                
                param.kerType = 2;
                param.rfeC = 2^n(acc_maxid(find(ismember(acc_maxid,MSE_minid),1)));
                param.useCBR = 1;
%                 % Feature Selection                
                [ftRank,ftScore] = ftSel_SVMRFECBR(double(training),double(trainlab(t,:))',param);
                ridx = find(sum(~isnan(ftScore),1) > (1/10)*size(ftScore,1));
%                 
                if 1
                    pixelmap = nan(xpix/p,ypix/p);
                    keeppix = 1:prod(size(pixelmap));
                    keeppix(throwpix) = [];
                    keeppix(rmpixel) = [];
                    pixelmap(keeppix)= nanmean(ftScore(:,ftRank),1);
                    subplot(2,nrtask,(t-1)*2+1); imagesc(pixelmap')
                    axis square
                    ylabel(['task '  taskOptNames{t}])
                    title(['Mouse ' mouse ' TW' num2str(TWidx) 'Fold ' num2str(nc) 'Task = ' num2str(t)])
                    pixelmap = nan(xpix/p,ypix/p);
                    valmat=length(ftRank):-1:1;
                    valmat=valmat(ftRank);
                    pixelmap(keeppix) = valmat;
                    subplot(2,nrtask,(t-1)*2+2); imagesc(pixelmap')
                    axis square
                    title('The higher the more informative')
                    
                    drawnow
                end
                
                %%  Actual prediction
                model = svmtrain(double(trainlab(t,:))',double(training(:,ftRank(ridx))),['-q -c ' num2str(2^n(acc_maxid(find(ismember(acc_maxid,MSE_minid),1))))]);
                [lbl, acc, dec] = svmpredict(double(testlab(t,:))',double(test(:,ftRank(ridx))),model,['-q']);
%                 model = svmtrain(double(trainlab(t,:))',double(training),['-q -c ' num2str(2^n(acc_maxid(find(ismember(acc_maxid,MSE_minid),1))))]);
%                 [lbl, acc, dec] = svmpredict(double(testlab(t,:))',double(test),model,['-q']);
%                 %
                %                 supVec = full(model.SVs);
                %
                %
                %                 alpha_signed = model.sv_coef;
%                 nSv = size(supVec,1);
%                 svInProd = supVec*supVec';
%                 svSqr = sum(supVec.^2,2);
%                 kerMatAll0 = repmat(svSqr,1,nSv) + repmat(svSqr',nSv,1) - 2*svInProd;
%                 w2_allIn = trace(alpha_signed' * exp(-2^-6 * kerMatAll0) * alpha_signed);
%                 % trace is used to add up the feature weights of each binary-class
%                 % subproblems. This strategy hasn't been verified.
%                 nFtIn = length(ridx);
%                 w2_in = zeros(1,nFtIn);
%                 
%                 for iFtIn = 1:nFtIn
%                     supVecP = supVec(:,iFtIn);
%                     % use the method in spider toolbox to compute the weight for
%                     % each feature p.
%                     kerMatP0 = (repmat(supVecP,1,nSv) - repmat(supVecP',nSv,1)).^2;
%                     kerMatRemoveP = exp(-2^-6 * (kerMatAll0 - kerMatP0));
%                     
%                     % the approximate margin when feature p is removed
%                     % it can be proved that w2=alpha_signed'*K*alpha_signed=
%                     % sum(alpha) when alpha is approximated
%                     w2_in(iFtIn) = trace(alpha_signed'*kerMatRemoveP*alpha_signed);
%                     % 				if rem(p,100)==0, fprintf(','); end
%                 end
%                 weightvec = nan(1,length(ftRank));
%                 weightvec(ftRank(ridx)) = w2_in; 
%                 figure; 
%                 pixelmap = nan(xpix/p,ypix/p);
%                 keeppix = 1:prod(size(pixelmap));
%                 keeppix(throwpix) = [];
%                 pixelmap(keeppix)= weightvec;
%                 h = imagesc(pixelmap');
%                 set(h,'alphadata',~isnan(pixelmap'))
                
                %                         flag = 0;
                %                         if length(unique(lbl)) == 1
                %                             warning('All labels the same..')
                %                             flag = 1;
                %                         end
                %                         if sum(double(nested_testl)==1)./length(double(nested_testl)) == acc(1) || sum(double(nested_testl)==2)./length(double(nested_testl))
                %                             disp('Exact same accuracy as number of labels...')
                %                             flag = 1;
                %                         end
                flag = 0;
                if length(unique(lbl)) == 1
                    warning('All labels the same..')
                    flag = 1;
                end
                if sum(double(testlab(t,:))==1)./length(double(testlab(t,:))) == acc(1)/100  ||  sum(double(testlab(t,:))==1)./length(double(testlab(t,:))) == (100-acc(1))/100
                    disp('Exact same accuracy as number of labels...')
                    flag = 1;
                end
                try
                    ACCURACY(t,nc) = acc(1);
                    MSE(t,nc) = acc(2);
                    SAMEWARNING(t,nc) = flag;
                    try
                    FeatureRank{t,nc} = ftRank;
                    FeatureScore{t,nc}= ftScore;
                    catch ME
                       disp(ME)
                    end
                    PARAMS{t,nc} = param;
                    MODELS{t,nc} = model;
                    
                catch ME
                    disp(ME)
                    keyboard
                end
                     
                actvspred(t,testcidx,:) = [testlab(t,:)',lbl];
                disp(['Task ' num2str(t) ', TW ' num2str(TWidx) ', fold ' num2str(nc) ' completed'])
                
            end
            
        catch ME
            disp(ME)
            keyboard
        end
        
      
    end    
    Results.ACCURACY = ACCURACY;
    Results.MSE = MSE;
    Results.SAMEWARNING = SAMEWARNING;
    Results.ActualVSPrediction = actvspred;
    try
    Results.FeatureRank = FeatureRank;
    Results.FeatureScore = FeatureScore;
    catch ME
        disp(ME)
    end
    Results.PARAMS = PARAMS;
    Results.ThrowPix = throwpix;
    Results.rmpixel = rmpixel;
    Results.MODELS = MODELS;
    %Save Results
    save(fullfile(ResultPath, [mouse 'time' num2str(TW{TWidx}(1)) '-' num2str(TW{TWidx}(2)), Model.Rnames{regio2take(roidx)}]),'Results')
    
    
end
end
