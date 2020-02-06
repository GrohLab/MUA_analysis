

clear all
files={'Z:\Ross\Experiments\17.12.19\Forepaw VPL\M167_L6-Stim10mW_Mech_CFA_ForePaw_Conc'
    'Z:\Ross\Experiments\13.11.19\M160_L6Stim+Mech_Saline'
    'Z:\Ross\Experiments\18.12.19\M168_10mW_MechStim_Saline_VPL'}

I=1;
[FILEPATH,NAME,EXT] = fileparts(files{I});

cd(FILEPATH)
load([NAME 'analysis.mat'],'Conditions','Triggers')
load([NAME '_sampling_frequency.mat'])
load([NAME '_CondSig.mat'],'chan67','head67')
%load(files{I})



%%
close all
clear params EEG spikeFindingData Conditions Triggers
load(files{I},'EEG','spikeFindingData','Conditions','Triggers')
LFP=chan67;
ppms=fs/1000;

%% all of this general enough for both old and new
%all triggers for spontaneous:
trigs=[];
if isstruct(Conditions)  %make compatible with E's more sensible Conditions variable
    mech=Conditions(1).Triggers;
    for j=1:numel(Conditions);tr=Conditions(j).Triggers(:,1); if size(tr,1)>size(tr,2), tr=tr';end;trigs=[trigs tr];end
else
    mech=Conditions{1}.Triggers;
    for j=1:numel(Conditions);tr=Conditions{j}.Triggers; if size(tr,1)>size(tr,2), tr=tr';end;trigs=[trigs tr];end
end

%stick signals together
cdata=double([LFP';Triggers.whisker';Triggers.laser']);
%gather spont epochs===================================
timeBefore=-5;timeAfter=-.5;
X=TriggeredEpochs(cdata,trigs,fs*timeBefore, fs*timeAfter);
lfp_spont=squeeze(X(1,:,:));
lfp_spont=lfp_spont(:,~isnan(sum(lfp_spont))); %get rid of bad trials)
%lfp_spont=resample(lfp_spont,1,round(ppms));%Downsample to 1kHz
lfp_spont=lfp_spont(1:round(ppms):end,:);%Downsample to 1kHz

%mechanical stimulus
timeBefore=.5;timeAfter=5;
X_mech=TriggeredEpochs(cdata,mech,fs*timeBefore, fs*timeAfter);
lfp_mech=squeeze(X_mech(1,:,:));
lfp_mech=lfp_mech(:,~isnan(sum(lfp_mech))); %get rid of bad trials)
lfp_mech=resample(lfp_mech,1,round(ppms));%Downsample to 1kHz
mechstim=resample(squeeze(X_mech(2,:,:)),1,round(ppms));


%%
figure
plot(mechstim)  %check the stimulus
%CHRONUX PARAMETERS
params.Fs=fs/ppms;  %WE DOWNSAMPLE af
params.fpass=[0 120]; % band of frequencies to be kept
Ktapers=5;
NW=(Ktapers+1)/2;
params.tapers=[NW Ktapers]; % taper parameters
params.pad=2; % pad factor for fft
params.err=[2 0.05];
params.trialave=1;
movingwin=[.5 0.02];
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