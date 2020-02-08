

%analyze LFP to compare CFA and Saline


clear all
files ={'Z:\Sailaja\Manuscript\Data\Juxtasomal\CFA\M107_C2\M107_C2_Mech+5mW_UMS_analysis.mat'...
    'Z:\Sailaja\Manuscript\Data\Juxtasomal\CFA\M110_C4\M110_C4_Mech+6mW_UMS_analysis.mat'...
    'Z:\Sailaja\Manuscript\Data\Juxtasomal\CFA\M110_C5\M110_C5_Mech+6mW_UMS_analysis.mat'...
    'Z:\Sailaja\Manuscript\Data\Juxtasomal\CFA\M113_C4\M113_C4_UMS_analysis.mat' ... %good?
    'Z:\Sailaja\Manuscript\Data\Juxtasomal\CFA\M113_C6\M113_C6_Mech+L6_UMS_analysis.mat'...
    'Z:\Sailaja\Manuscript\Data\Juxtasomal\CFA\M133_C2\M113_C2_Mech+L6alone edited_UMS_analysis.mat'...
    'Z:\Sailaja\Manuscript\Data\16 channel\Poly2A probe\CFA\L6 alone\27_12_2018\M47_C2_L6Mech+L61mWanalysis.mat'...
    %    'Z:\Sailaja\Manuscript\Data\16 channel\Poly2A probe\27_12_2018\M47_C2_L6Mech+L61mWanalysis.mat'...  %this is a repeat?

    Z:\Sailaja\Manuscript\Data\16 channel\Poly2A probe\CFA\VPL and L6\09_11_2018\M3_C1   M3_C1_L6_Mechanalysis.mat
  %  Z:\Sailaja\Manuscript\Data\16 channel\Poly2A probe\CFA\VPL and L6\09_11_2018\M3_C2  %no analysis!
  
  Z:\Sailaja\Manuscript\Data\16 channel\Poly2A probe\CFA\VPL and L6\12_3_2019\M59_C1  M59_C1_HL +Terminal sti_1mWanalysis.mat
  Z:\Sailaja\Manuscript\Data\16 channel\Poly2A probe\CFA\VPL and L6\12_3_2019\M59_C2\C2A M59_C2_HL +Terminal sti_1mWanalysis.mat
  Z:\Sailaja\Manuscript\Data\16 channel\Poly2A probe\CFA\VPL and L6\12_3_2019\M59_C2\C2B M59_C2_HL +Terminal sti_2mWanalysis.mat
  
    'Z:\Sailaja\Manuscript\Data\16 channel\Poly2A probe\M7_C3\M7_C3_Mech_05mWanalysis.mat'...
    
    
    'Z:\Sailaja\Manuscript\Data\16 channel\standard 16 channel\probe\M137_C4\M137_C4_Mech+L6 05mWanalysis.mat'}

%% ROSS FILES HERE

% Saline VPL
    % 'D:\Ross\18.12.19' - 10mW L6 With LFP             Counts by condition
        % on SDS
    % 'D:\Ross\13.11.19' - 3.5mW L6 (Don't think LFP)   Counts by condition
        % on SDS
    % 'D:\Ross\4.02.20' 2.5mW - Kilosorting but should be curated by end of 08.02.20
                
         
% CFA VPL
    % 'D:\Ross\17.12.19\Forepaw VPL' With LFP       Counts by condition
        % on SDS
    % 'D:\Ross\30.01.20' With LFP
    
% CFA Cortex
    
    % 'D:\Ross\17.12.19\Hind Limb Cortex'   Counts by condition
        % on SDS
        
% DREADD Expts

    % Gi 'D:\Ross\10.12.19\Conc'
    % Gq 'D:\Ross\11.12.19\Concatenated'

%%

for I=9:numel(files)
    close all
    clear params EEG spikeFindingData Conditions Triggers
    load(files{I},'EEG','spikeFindingData','Conditions','Triggers')
    lfp=EEG.data;
    ppms=spikeFindingData.ppms;
    fs=spikeFindingData.ppms*1000;
    mech=Conditions{1}.Triggers;
    %all triggers for spontaneous:
    trigs=[];
    for j=1:numel(Conditions);tr=Conditions{j}.Triggers; if size(tr,1)>size(tr,2), tr=tr';end;trigs=[trigs tr];end
    
    %lfp=zscore(EEG.data);
    
    
%     %simple spectrogram: spontaneous
%     %chronux parameters;
%     Params.Fs=fs  %WE WILL DOWNSAMPLE
%     params.fpass=[0 200]; % band of frequencies to be kept
%     Ktapers=12;
%     NW=(Ktapers+1)/2;
%     params.tapers=[NW Ktapers]; % taper parameters
%     params.pad=2; % pad factor for fft   

    % get rid line noise (50 Hz)
    %[f] = delineTest(lfp,fs)
   % LFP = delineSignal(lfp,fs,(1:15)*f, 4);  %what is 1:9??? and 4???
    LFP=lfp;
    
     
    %stick signals together
    cdata=[LFP';Triggers.whisker';Triggers.light'];
    
   
    
    %gather stim-triggered epochs===================================
    timeBefore=-10.5;timeAfter=-.5;
    X=TriggeredEpochs(cdata,trigs,fs*timeBefore, fs*timeAfter);
    lfp_spont=squeeze(X(1,:,:));
    lfp_spont=lfp_spont(:,~isnan(sum(lfp_spont))); %get rid of bad trials)
    lfp_spont=resample(lfp_spont,1,round(ppms));%Downsample to 1kHz
    
    %mechanical stimulus
    timeBefore=-1;timeAfter=5;
    X_mech=TriggeredEpochs(cdata,mech,fs*timeBefore, fs*timeAfter);
    lfp_mech=squeeze(X_mech(1,:,:));
    lfp_mech=lfp_mech(:,~isnan(sum(lfp_mech))); %get rid of bad trials)
    lfp_mech=resample(lfp_mech,1,round(ppms));%Downsample to 1kHz
    mechstim=resample(squeeze(X_mech(2,:,:)),1,round(ppms));
    
    figure
    plot(mechstim)  %check the stimulus
    
    
  
    %CHRONUX PARAMETERS
    params.Fs=fs/ppms;  %WE DOWNSAMPLE af
    params.fpass=[0 120]; % band of frequencies to be kept
    Ktapers=12;
    NW=(Ktapers+1)/2;
    params.tapers=[NW Ktapers]; % taper parameters
    params.pad=2; % pad factor for fft
    params.err=[2 0.05];
    params.trialave=1;
    movingwin=[2 0.05];
    segave=1;
    winseg=1
    
    
    % power spectra for spontaneous and mechanical
    lfp_spont=rmlinesc(lfp_spont,params); %remove line noise for ds signal
    lfp_mech=rmlinesc(lfp_mech,params); %remove line noise for ds signal
    
    figure   %plot it all
    [S,f]=mtspectrumc(lfp_spont,params);
    plot(f,10*log10(S)); xlabel('Frequency Hz'); ylabel('Spectrum');
    spectra.Sspont=S; spectra.f=f;
    hold on
    [S,f]=mtspectrumc(lfp_mech,params);
    plot(f,10*log10(S)); xlabel('Frequency Hz'); ylabel('Spectrum');
    spectra.Smech=S;

    %     figure;
    %     [S,f, Serr]=mtspectrumc(lfp_trials,params);
    %     plot(f,10*log10(S),f,10*log10(Serr(1,:)),f,10*log10(Serr(2,:))); xlabel('Frequency Hz'); ylabel('Spectrum');
    %%% pause
    
    
    %time resolved spectrogram for before and after stimulus trigger (t=0)
    timeBefore=-5;timeAfter=12; 
    condSpect=struct([]) 
    for ii=1:numel(Conditions)
        trigs=Conditions{ii}.Triggers;
        trigs=trigs(find(trigs<(size(cdata,2)-timeAfter*fs)));
        X=TriggeredEpochs(cdata,trigs,fs*timeBefore, fs*timeAfter);
        lfp_trials=squeeze(X(1,:,:));
        lfp_trials=lfp_trials(:,~isnan(sum(lfp_trials))); %get rid of bad trials)
        lfp_trials=resample(lfp_trials,1,round(ppms));%Downsample to 1kHz
        lfp_trials=rmlinesc(lfp_trials,params);%remove line noise for ds signal

        stim=squeeze(X(2,:,:));
        stim=stim(:,~isnan(sum(stim))); %get rid of bad trials)
        stim=resample(stim,1,round(ppms));%Downsample to 1kHz
        tt=[(timeBefore*fs):(timeAfter*fs)]/fs;
        tt=resample(tt,1,round(ppms));%Downsample to 1kHz

    
        light=squeeze(X(3,:,:));
        light=light(:,~isnan(sum(light))); %get rid of bad trials)
        light=resample(light,1,round(ppms));%Downsample to 1kHz
        tt=[(timeBefore*fs):(timeAfter*fs)]/fs;
        tt=resample(tt,1,round(ppms));%Downsample to 1kHz
        
        [S,t,f]=mtspecgramc(lfp_trials,movingwin,params);
        S=S';
        condSpect(ii).tt=tt;
        condSpect(ii).t=t;
        condSpect(ii).f=f;
        condSpect(ii).params=params;
        condSpect(ii).S=S;
        condSpect(ii).stim=mean(stim'/max(max(stim)));
        condSpect(ii).light=mean(light'/max(max(light)))
    end
    

    blwindow=(t<=-timeBefore);
    baseline=mean(S(:,blwindow)');
    dbS=10*log10((S./baseline'));
    figure
    pcolor(t+timeBefore,f',dbS);shading flat;
    %pcolor(t+timeBefore,f',S);shading flat;

    spectrogram.S=S; spectrogram.f=f;spectrogram.t=t;
    spectrogram.params=params;
  

    colormap fire
    hold on
    plot(tt,mean(stim'/max(max(stim)))*10 - 11);
    ylim([-10 120])
    %plot(zeros(size(f)),f,'-k','linewidth',2)
    [FILEPATH,NAME,EXT] = fileparts(files{I}); 
    cd(FILEPATH)
    save ChronuxResults 'spectrogram' 'spectra' Conditions condSpect NAME FILEPATH

    
end
  
%% SCRAP
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


[imf,residual,info] = emd(LFP,'Interpolation','pchip', 'MaxNumIMF',20);
[hs,f,t,imfinsf,imfinse]=hht(imf,fs/round(ppms),'FrequencyLimits',[0 200]);
goods=find(max(imfinsf)>3);
new_lfp=sum(imf(:,goods)');
figure
plot(LFP)
hold on;plot(new_lfp)

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




