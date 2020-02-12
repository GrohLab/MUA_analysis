 %% multiunit recording practice
clearvars
% cd 'D:\Dropbox\16 Channel Recording may 2018'
% homedir='F:\Experiments_2018\16 channel\Standard probe\19_4_2018\M137_C5';
homedir='D:\Ross\18.12.19';
%fname = '17.12.19';
expName = 'M168_10mW_MechStim_Saline_VPL';
fname = expName;
load(fullfile(homedir,[fname,'_all_channels.mat']))
load(fullfile(homedir,[fname,'analysis.mat']),'Conditions','Triggers')
if ~iscell(Conditions)
    auxCell = cell(numel(Conditions),1);
    for ccond = 1:numel(Conditions)
        auxCell(ccond) = {Conditions(ccond)};
    end
    Conditions = auxCell;
    clearvars auxCell
end
Ncond = numel(Conditions);
if exist('Fs','var')
    fs = Fs;
elseif exist('fs','var')
    Fs = fs;
elseif exist([fname,'_sampling_frequency.mat'],'file')
    try
        load([fname,'_sampling_frequency.mat'],'fs')
    catch
        % nothing happens
    end
else
    strFS = inputdlg('Unknown sampling frequency. Please provide it:'...
        ,'Sampling Frequency warning');
    fs = str2double(strFS);
    save(fullfile(homedir,[fname,'_sampling_frequency.mat']),'fs')
end
ppms=fs/1000;


%% Initialize the variables
% This section loads the cluster spike times into the 'Spikes' cell array.
% If there is a step gone wrong, you can re-initialize re-running this
% section.

%population PSTHs
Ntcl = size(sortedData,1);
Spikes=cell(Ntcl,1);
Names=cell(Ntcl,1);

mech = Triggers.whisker;
mObj = StepWaveform(mech,fs);
mSubs = mObj.subTriggers;
mech = mObj.subs2idx(mSubs,mObj.NSamples);
try
    light = Triggers.light;
catch
    light = Triggers.laser;
end
lObj = StepWaveform(light,fs);
lSubs = lObj.subTriggers;
light = lObj.subs2idx(lSubs,lObj.NSamples);
continuousSignals = {mech;light};
Ns = length(mech);
Nt = Ns/fs;

for i=1:size(sortedData,1)
    if iscolumn(sortedData{i,2})
        Spikes{i}=cell2mat(sortedData(i,2))';
    else
        Spikes{i}=cell2mat(sortedData(i,2));
    end
    Names{i}=sortedData(i,1);
end
badsIdx = cellfun(@(x) x==3,sortedData(:,3));
bads = find(badsIdx);
totSpkCount = cellfun(@numel,sortedData(:,2));
silentUnits = totSpkCount/Nt < 0.1;
bads = union(bads,find(silentUnits));
goods=setdiff(1:size(sortedData,1),bads);

badsIdx = StepWaveform.subs2idx(bads,size(sortedData,1));

%% DATA EXPLORER SECTION
% Viewing window and bin size both in seconds
timeLapse = [1, 6];
binSz = 0.02;
Ngc = numel(goods);
% Logical trace for the first considered cluster and column sample
% subscripts for the rest.
spkLog = StepWaveform.subs2idx(round(sortedData{goods(1),2}*fs),Ns);
spkSubs = cellfun(@(x) round(x.*fs),sortedData(goods(2:end),2),'UniformOutput',false);
fig = gobjects(Ncond,1);
for conSel = 1:Ncond
    if contains(Conditions{conSel}.name,'all','IgnoreCase',true)
        continue
    end
    [dst, cst] = getStacks(spkLog,Conditions{conSel}.Triggers,'on',timeLapse,fs,...
        fs,spkSubs,continuousSignals);
    stims = mean(cst,3);
    [Ne, Nt, Na] = size(dst);
    [PSTH, trig, sweeps] = getPSTH(dst,timeLapse,false(Na,1),binSz,fs);
    IDe = [Conditions{conSel}.name;sortedData(goods,1)];
    fig(conSel) = plotClusterReactivity(PSTH,trig,sweeps,timeLapse,binSz,...
        IDe,sprintf('%s\n%s',...
        expName, Conditions{conSel}.name),stims,{'Piezo','Laser'});
    configureFigureToPDF(fig(conSel));
    print(fig(conSel),fullfile(homedir,sprintf('%s %s.pdf',...
        expName, Conditions{conSel}.name)),'-dpdf','-fillpage')
end


%% look at inter-spike intervals as another check for pure light-evoke electrical artifacts
close all
% for i=1:numel(Spikes) % insert cluster numbers here
%     figure
%     isis=[diff(Spikes{i})*1000]; %in ms
%     [hisi isi_bins]=hist(isis,[1:10000]);
%     plot(log10(isi_bins),hisi,'LineWidth',2)
%     title(['Cluster: ',sortedData{i,1},' #',num2str(i)])
% end
spkSubs = cellfun(@(x) round(x.*fs),sortedData(goods,2),'UniformOutput',false);
[~,FRpu] = plotISI(spkSubs,fs,sortedData(goods,1));


%%

%bads=[1 2 3 14 6]  %1 2 3 14 are light artifacts
%bads=[2, 3, 4, 5, 6, 7, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 21,...
%    22, 23, 25, 26, 27, 31];

noresponse=[];
bads=[bads noresponse];
goods = setdiff(1:size(sortedData,1),bads);
%bads=[]; %uncomment this to have bads empty
%assumes that all spike trains are good


%removes bad units


%Spikes=Spikes(goods);
%Names=Names(goods);




%% looking at collected data

close all
timeBefore=1000*ppms;
timeAfter=6000*ppms;
binsize=2*ppms;
plotit=1;
%spikes back in samples
try
    spikes=cell2mat(Spikes).*1000.*ppms; 
catch
    spikes=cell2mat(Spikes').*1000.*ppms; 
end
for I=1:numel(Conditions)
    name=Names{I};
    if size(Conditions{I}.Triggers,2) > 1
        triggers=Conditions{I}.Triggers(:,1);
    else
        triggers=Conditions{I}.Triggers;
    end
    [sp h bins trig_mech trig_light f]=triggeredAnalysisMUA(spikes,ppms,triggers,binsize,timeBefore, timeAfter,Conditions{I}.name,mech,light,plotit);
    fig = gcf;
    ax = fig.Children;
    linkaxes(ax,'x');
    title(Conditions{I}.name)
end
%% Cross-correlation or relationship between clusters
% Creating a square matrix containing the cross correlation normalized
% peaks omitting the 'bads' clusters. The input to the cross correlation
% function should be 'normalized'. This means that the spike traces shall
% be reconstructed from the time stamps and input to the xcorr function.
lenSpks = length(Spikes);
consIdxs = true(1,lenSpks);
% Some clusters were eliminated after an induvidual inspection of their
% PSTHs.
% bads = [];
consIdxs(bads) = false;
crscor = zeros(lenSpks,lenSpks,3);    % Square matrix with Ncl x Ncl
for ccl = 1:lenSpks
    xccl = ccl+1;
    % Cross correlation avoiding the autocorrelations and the 'bads'.
    while consIdxs(ccl) && xccl <= lenSpks
        % auxXCOR = xcorr(Spikes{ccl},Spikes{xccl},'coeff');
        if consIdxs(xccl)
            clstrInfo =...
                ['Cluster ',num2str(ccl),' against cluster ',num2str(xccl)];
            disp(clstrInfo)
            % Cross correlation
            % [auxCorr, lTx] = xcorr(auxSignal1,auxSignal2,'coeff');
            % figure;plot(lTx/fs,auxCorr)
            % Distance matrix
            dfMtx = log(distmatrix(Spikes{ccl}',Spikes{xccl}')+1);
            lnIdx = dfMtx < log(16);
            [y,x]=find(lnIdx);
            [mdl,yhat,r2] = fit_poly(x,y,1);
            eqInfo = ['y = ',num2str(mdl(1)),'x ',num2str(mdl(2))];
            display(eqInfo)
            figure('Name',clstrInfo);
            imagesc(dfMtx);hold on;plot(x,yhat,'LineStyle','--',...
                'LineWidth',3,'Color',[1,0,0]);title([eqInfo,' ',num2str(r2)])
            crscor(ccl,xccl,1:2) = mdl;
            crscor(ccl,xccl,3) = r2;
        end
        xccl = xccl + 1;
    end
end

%save cross-correlation so that we don't need to redo this constantly
cd(homedir)
save(fullfile(homedir,'CrossCoeffData.mat'),'crscor','goods','bads')

%% Merge similar clusters
% The marging packages indicate which clusters should be merged together
% due to their high similarity. The possibilities are that they belong to a
% same unit as busrting spikes, the cell shifted to another channel or any
% other reasonable cause.
mergingPackages = {[]};
Npg = numel(mergingPackages);
mSpikes = cell(1,Npg);
auxSignal = false(1,length(mech));
for cpg = 1:Npg
    for ccl = 1:numel(mergingPackages{cpg})
        auxSignal(round(fs*Spikes{mergingPackages{cpg}(ccl)})) = true;
    end
    mSpikes(cpg) = {find(auxSignal)/fs};
    auxSignal = false(1,length(mech));
    bads = [bads mergingPackages{cpg}(2:end)];
    Spikes(mergingPackages{cpg}(1)) = {single(mSpikes{cpg})};
end
bads = sort(unique(bads));

%save merging info
save CrossCoeffData Spikes mergingPackages bads -append


%% looking at individual conditions and clusters for good clusters

plotRasterFlag = true; % No raster plot in the figures!!
plotit = true;
timeBefore=250*ppms;
timeAfter=350*ppms;
binsize=2*ppms;
%close all
compTime = 0.11;
preTrigWindow = [-compTime,0] * 1e3;
posTrigWindow = [0, compTime] * 1e3;
inactivUnits = zeros(size(goods,1),1,'single');
respUnits = zeros(size(goods,1),1,'single');
cgu = 0;
for i=goods
    if ~ismember(i,bads)
        for I=1:4 %pick out conditions to look at
            cgu = cgu + 1;
            spikes=(Spikes{i})*1000*ppms; %spikes back in samples
            name=Names{i};
            if size(Conditions{I}.Triggers,2) > 1
                triggers=Conditions{I}.Triggers(:,1);
            else
                triggers=Conditions{I}.Triggers;
            end
            Na = length(triggers);
            [sp h bins trig_mech trig_light f]=...
                triggeredAnalysisMUA(spikes,ppms,triggers,binsize,timeBefore, timeAfter,Conditions{I}.name,mech,light,plotit,plotRasterFlag);
            inactivUnits(cgu) = sum(cellfun(@isempty,sp))/numel(triggers);
            
            Nts = sum(cellfun(@numel,sp));
            spikeDensity = Nts/Na;
            fig = gcf;
            if inactivUnits(cgu) > 0.7 || spikeDensity < 0.3
                close(fig)
                % Cluster number, cluster ID, spike density, inactive
                % trials
                fprintf(1,'%d ID:%s SD:%.2f AT:%.2f%%\n',i,sortedData{i,1},...
                    spikeDensity,(1-inactivUnits(cgu))*100)
                continue
            end
            h(h==0) = min(h(h~=0))*1e-6;
            preTrig = h(bins > preTrigWindow(1) & bins < preTrigWindow(2))*...
                Na;
            posTrig = h(bins >= posTrigWindow(1) & bins < posTrigWindow(2))*...
                Na;
            respUnits(i) = sum(preTrig)/numel(preTrig) -...
                sum(posTrig)/numel(posTrig);
            title(sprintf('Cluster: %s #%d A: %.2f E: %.3f R: %.3f',...
                sortedData{i,1},i,(1-inactivUnits(cgu))*100,getEntropyFromPDF(h),...
                respUnits(i)))
            ax = fig.Children;
            ax(2).Title.String = sprintf('%.2f av. spikes per trial',spikeDensity);
            linkaxes(ax,'x');            
        end
        
    end
end


%% looking at possible artifacts
wru = [1, 6, 13, 15, 16, 17, 21, 38, 40, 43, 45];
timeBefore=250*ppms;
timeAfter=350*ppms;
binsize=2*ppms;
plotit = true;
plotRasterFlag = true; % raster plot in the figures!!
% 16 Gaussians, 3 parameters estimated, 7 conditions, 8 clusters
Ngauss = 16;
Ncond = 7;
Nclus = length(wru);
paramsCond = nan(Ngauss,3,Ncond,Nclus);
cclu = 1;
responseWindow = [0, 0.25]; % In seconds
xdom = responseWindow(1):1/fs:responseWindow(2);
E = 1e-7;
for i= wru % insert cluster numbers here
    for I = 3:9
        spikes=(Spikes{i})*1000*ppms; %spikes back in samples
        name=Names{i};
        if size(Conditions{I}.Triggers,2) > 1
            triggers=Conditions{I}.Triggers(:,1);
        else
            triggers=Conditions{I}.Triggers;
        end
        [sp h bins trig_mech trig_light f]=...
            triggeredAnalysisMUA(spikes,ppms,triggers,binsize,...
            timeBefore,timeAfter,Conditions{I}.name,mech,light,...
            plotit,plotRasterFlag);
        spJoint = cell2mat(sp')/fs;
        spJoint = spJoint(spJoint >= responseWindow(1) &...
            spJoint <= responseWindow(2));
        if ~isempty(spJoint) && numel(spJoint) > 2
            paramsAux = emforgmm(spJoint, Ngauss, E, false);
            M = size(paramsAux,1);
            paramsCond(1:M,:,I-2,cclu) = paramsAux;
        else
            fprintf(1,'No response in cluster %d, condition %s',...
                i,Conditions{I}.name)
            fprintf(1,' during %.3f and %.3f ms\n',responseWindow(1)*1e3,...
                responseWindow(2)*1e3);
        end
        if plotit
            fig = gcf;
            ax = fig.Children;
            linkaxes(ax,'x');
            
            title(['Clus: ',num2str(i),' Cond: ',Conditions{I}.name])
        end
    end
    cclu = cclu + 1;
end
%% Spikes backup.
% Elimination of the 'bads' in the 'Spikes' variable.
% But also making a backup before deleting.
% WARNING! If you run this section twice without restoring the backup, you
% will lose good units!!
Spikes_BACKUP = Spikes;
Spikes(bads) = [];
consIdxs(bads) = [];
%% Restore Spikes
% Use this section to restore the 'bads' labelled units into the 'Spikes'
% variable
Spikes = Spikes_BACKUP;
%% organize real figure


% population histogram for later plotting, by condition
spikes=cell2mat(Spikes)*1000*ppms;
plotit=true;
timeBefore=1000*ppms;
timeAfter=6000*ppms;
binsize=2*ppms;
H=[];
conds={};
count=0;
Trig_mech={};Trig_light={};
Ngauss = 16;
Ncond = 7;
Nclus = length(wru);
paramsCond = zeros(16,3,2);

cc = 1;
responseWindow = [0, 0.1]; % In seconds
for I=3:numel(Conditions)
    count=count+1;
    %spikes back in samples
    name=Names{I};
    if size(Conditions{I}.Triggers,2) > 1
        triggers=Conditions{I}.Triggers(:,1);
    else
        triggers=Conditions{I}.Triggers;
    end
    conds{count}=Conditions{I}.name;
    [sp, h, bins, Trig_mech{count}, Trig_light{count}, f]=...
        triggeredAnalysisMUA(spikes,ppms,triggers,binsize,timeBefore,...
        timeAfter,Conditions{I}.name,mech,light,plotit);
    spJoint = cell2mat(sp')/fs;
    spJoint = spJoint(spJoint >= responseWindow(1) &...
        spJoint < responseWindow(2));
    paramsCond(:,:,cc) = emforgmm(spJoint,16,1e-10,0);
    cc = cc + 1;
    title(Conditions{I}.name)
    fig = gcf;
    ax = fig.Children;
    linkaxes(ax,'x')
    %convert to rate in Hz
    h=h*(1000/binsize*ppms);
    H(count,:)=h;
end
t=[-timeBefore:timeAfter]/ppms;
%
%% get individual responses by condition

Sp={};
% close all
count=0;
plotit=0;

for I=1:numel(Conditions)  %by condition
    count=count+1;
    if size(Conditions{I}.Triggers,2) > 1
        triggers=Conditions{I}.Triggers(:,1);
    else
        triggers=Conditions{I}.Triggers;
    end
    SPIKES={};
    for i=1:numel(Spikes)   %for each neuron
        if consIdxs(i)
            spikes=(Spikes{i})*1000*ppms; %spikes back in samples
            name=Names{i};
            %use same binning etc as pop hist
            [SPIKES{i} h bins trig_mech trig_light f]=triggeredAnalysisMUA(spikes,ppms,triggers,binsize,timeBefore, timeAfter,Conditions{I}.name,mech,light,plotit);
        end
    end
    Sp{count}=SPIKES;
end



%

%% get data appropriate for plotting rasters


SPIKESs={};YSs={}; %all data here ,cond x num neuron

for ii=1:numel(Sp) %pick one condition
    sp=Sp{ii};
    shift=0; SPIKES={};YS={};
    for j=1:numel(sp) %over all neurons
        spikes=sp{j};
        ys={};
        for jj=1:numel(spikes) %over all trials
            if ~isempty(spikes{jj})
                ys{jj}=ones(size(spikes{jj}))+shift;
            else
                ys{jj}=[];
            end
            shift=shift+1;%add for each trial, per neuron
        end
        SPIKES{j}=cell2mat(spikes');
        YS{j}=cell2mat(ys');
        
    end
    
    SPIKESs{ii}=SPIKES;
    YSs{ii}=YS;
end
%

%% get color for each neuron, just for plotting
% cmap=jet();
% n=floor(size(cmap,1)/numel(SPIKESs{1}));
% colors=cmap(1:n:end,:);
colors = jet(numel(SPIKESs{1}));

%% plot it all
figure('Color',[1,1,1])
% yUpLimit = max(cell2mat(cellfun(@(x) (cellfun(@max,x(end),'UniformOutput',false)),YSs)));
Ncl = numel(SPIKES);
Ncon = numel(Conditions);
for ii=1:Ncon
    auxAx = subplot(6,Ncon,[ii, ii+Ncon, ii+(Ncon*2)]);
    Nt = numel(Conditions{ii}.Triggers);
    for j=1:numel(SPIKESs{ii})        
        xs=SPIKESs{ii}{j}/ppms;
        ys=YSs{ii}{j};
        plot(xs,ys,'.','color',colors(j,:),'MarkerSize',10)
        hold on
    end
    lstResp = find(~cellfun(@isempty,YSs{3}),1,'last');
    yUpLimit = max(YSs{ii}{lstResp});
    yticks = ((1:Ncl) - 0.5) * Nt;
    set(auxAx,'YTick',yticks,'YTickLabel',1:Ncl);ylabel(...
        sprintf('Clusters_{(t = %d)}',Nt))
    
    % ylabel 'trials/neuron'
    box off
    xlim([min(bins) max(bins)])
    ylim([0,yUpLimit])
    
end



% titles={'mechanical',...
%     'mechanical + 1 Hz L6',...
%     'mechanical + 10 Hz L6',...
%     '10 Hz L6 control'};
titles = cell(Ncon,1);
for cct = 1:Ncon
    titles(cct) = {Conditions{cct}.name};
end
yLimit = 5*ceil(5\(max(H(:)) * 1.05));
% yLimit = 50;
% load(fname,'chan21')
chan21 = reshape(mech,1,length(mech));
[~,cStack] =...
    getStacks(false(1,length(chan21)),...
    Conditions{1}.Triggers,'on',[-min(bins),max(bins)]*1e-3,...
    fs,fs,[],chan21);
meanMech = mean(squeeze(cStack),2);
for ii=1:Ncon
    
    auxAx = subplot(6,Ncon,[ii+(4*Ncon) ii+(5*Ncon)]);
    bar(bins,H(ii,:),'k')
    if ~strcmpi(Conditions{ii}.name,'lasercontrol')
        
    end
    xlim([min(bins) max(bins)])
    xlabel ms
    % ylabel 'pooled spike rate'
    title(titles{ii})
    box off
    ylim([0 yLimit])
end

%%
t = 0:1/fs:(length(Trig_mech{1}(1,:))-1)/fs;
t = 1000*(t - timeBefore/fs);
subplot(6,Ncon,4*Ncon - 3)
plot(t,Trig_mech{1}(1,:),'Color',[255, 128, 0]/255,'linewidth',2)
hold on;plot(t,...
    scaleOnLine(meanMech(1:length(t)),0,1),'Color',[255, 51, 0]/255,'linewidth',2)
set(gca,'Visible','off')
box off
xlim([min(bins) max(bins)])
ylim([0 1.5])
ylabel Stimulus

subplot(6,Ncon,4*Ncon - 2)
plot(t,Trig_mech{2}(1,:),'Color',[255, 128, 0]/255,'linewidth',2);
hold on;plot(t,...
    scaleOnLine(meanMech(1:length(t)),0,1),'Color',[255, 51, 0]/255,'linewidth',2)
hold on
plot(t,Trig_light{2}(1,:),'Color', [0, 64, 255]/255','linewidth',1);
set(gca,'Visible','off')
box off
xlim([min(bins) max(bins)])
ylim([0 1.5])

subplot(6,Ncon,4*Ncon - 1)
plot(t,Trig_mech{3}(1,:),'Color',[255, 128, 0]/255,'linewidth',2);
hold on;plot(t,...
    scaleOnLine(meanMech(1:length(t)),0,1),'Color',[255, 51, 0]/255,'linewidth',2)
hold on
plot(t,Trig_light{3}(1,:),'Color', [0, 64, 255]/255','linewidth',1);
set(gca,'Visible','off')
box off
xlim([min(bins) max(bins)])
ylim([0 1.5])

subplot(6,Ncon,4*Ncon)
plot(t,Trig_light{4}(1,:),'Color', [0, 64, 255]/255','linewidth',1);
set(gca,'Visible','off')
box off
xlim([min(bins) max(bins)])
ylim([0 1.5])

%% plot some raw data
figure
load 'M137_C5_Mech_L6 05mWanalysis','filteredResponse'

v=filteredResponse.data;indices=[10000:(1000*ppms*700)]+30*1000*ppms; v=v/max(v);
time=indices/ppms/1000;
light_in=light(indices);light_in=light_in/max(light_in);
mech_in=mech(indices);mech_in=mech_in/max(mech_in);
plot(time,v(indices)/2,'k')
hold on
for i=1:numel(Spikes)
    
    sp=Spikes{i}*1000*ppms;
    sp=sp(ismember(sp,indices));
    plot(sp/ppms/1000, ones(size(sp))+i*.05-.7,'.','color',colors(i,:),'markersize',15)
    hold on
    axis tight
    ylim([-1.4 1])
    
end
plot(time,light_in/4-.8,'c','linewidth',.5)
plot(time,mech_in/4-1.26,'g','linewidth',1.5)



