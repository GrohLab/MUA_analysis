

%analyze LFP to compare CFA and Saline


clear all
files ={'Z:\Sailaja\Manuscript\Data\Juxtasomal\CFA\M107_C2\M107_C2_Mech+5mW_UMS_analysis.mat'...
    'Z:\Sailaja\Manuscript\Data\Juxtasomal\CFA\M110_C4\M110_C4_Mech+6mW_UMS_analysis.mat'...
    'Z:\Sailaja\Manuscript\Data\Juxtasomal\CFA\M110_C5\M110_C5_Mech+6mW_UMS_analysis.mat'...
    'Z:\Sailaja\Manuscript\Data\Juxtasomal\CFA\M113_C4\M113_C4_UMS_analysis.mat' ... %good?
    'Z:\Sailaja\Manuscript\Data\Juxtasomal\CFA\M113_C6\M113_C6_Mech+L6_UMS_analysis.mat'...
    'Z:\Sailaja\Manuscript\Data\Juxtasomal\CFA\M133_C2\M113_C2_Mech+L6alone edited_UMS_analysis.mat'...
    'Z:\Sailaja\Manuscript\Data\16 channel\Poly2A probe\CFA\L6 alone\27_12_2018\M47_C2_L6Mech+L61mWanalysis.mat'...
    'Z:\Sailaja\Manuscript\Data\16 channel\Poly2A probe\M7_C3\M7_C3_Mech_05mWanalysis.mat'...
    'Z:\Sailaja\Manuscript\Data\16 channel\Poly2A probe\27_12_2018\M47_C2_L6Mech+L61mWanalysis.mat'...
    'Z:\Sailaja\Manuscript\Data\16 channel\standard 16 channel\probe\M137_C4\M137_C4_Mech+L6 05mWanalysis.mat'}

%% where is LFP for this?
%load('Z:\Sailaja\Manuscript\Data\16 channel\standard 16 channel probe\M137_C5\M137_C5_Mech_L6 05mWanalysis.mat','EEG','spikeFindingData','Conditions','Triggers')
%load('Z:\Sailaja\Manuscript\Data\16 channel\standard 16 channel probe\M137_C5\Old analysis_alex wrapper\M137_C5_Mech_L6 05mWanalysis.mat','EEG')


%% ROSS FILES HERE



%%
for I=2 %:numel(files)
    load(files{I},'EEG','spikeFindingData','Conditions','Triggers')
    ppms=spikeFindingData.ppms;
    fs=spikeFindingData.ppms*1000;
    mech=Conditions{1}.Triggers;
    lfp=zscore(EEG.data);
    
    cdata=[lfp';Triggers.whisker'/max(Triggers.whisker)];
    
    timeBefore=-3.5;timeAfter=10;
    X=TriggeredEpochs(cdata,mech,fs*timeBefore, fs*timeAfter);
    lfp_trials=squeeze(X(1,:,:));
    lfp_trials=double(resample(lfp_trials,1,round(ppms)));
    mechstim=squeeze(X(2,:,:));
    mechstim=double(resample(mechstim,1,round(ppms)));
    
    
    params.Fs=fs/round(ppms);
    params.fpass=[0 120]; % band of frequencies to be kept
    Ktapers=12;
    NW=(Ktapers+1)/2;
    params.tapers=[NW Ktapers]; % taper parameters
    
    params.pad=2; % pad factor for fft
    params.err=[2 0.05];
    params.err=0;
    params.trialave=1;
    movingwin=[2 0.05];
    segave=1;
    winseg=1
    
    [S,t,f]=mtspecgramc(lfp_trials,movingwin,params);
    
    S=S';
    blwindow=(t<=-timeBefore);
    baseline=mean(S(:,blwindow)');
    dbS=10*log10((S./baseline'));
    
    figure
    pcolor(t+timeBefore,f',dbS);shading flat;
    colormap fire
    hold on
    lfpdata.S=S;lfpdata.f
    plot(zeros(size(f)),f,'k-','linewidth',5)
   % save(files{I},'lfpdata','-append')
    clear lfpdata;
    
end
    
    %%
    LFP=resample(double(lfp),1,round(ppms)); LFP=LFP(1:400000);
    [imf,residual,info] = emd(LFP,'Interpolation','pchip');
    [hs,f,t,imfinsf,imfinse]=hht(imf,fs/round(ppms),'FrequencyLimits',[0 200]);
    
    %% HHT analysis  does not work for short recordings?
    figure
    ax=[];
    ax(1)=subplot(3,1,2)
    hht(imf,fs,'FrequencyLimits',[0 24]);
    colormap fire
    cmap=colormap;
    colormap(flipud(cmap))
    %%colormap(flipud(cb))
    %brighten(.5)
    %set(gca,'color',cmap(1,:))
    set(gca,'color',[1 1 1]*.7)
    
    %%
    
    ax(2)=subplot(3,1,3)
    % hht(imf,fs,'FrequencyLimits',[25 150],  'FrequencyResolution',1);
    % colormap fire
    % set(gca,'color',cmap(1,:))
   

