%matlabpoolv
clear all;close all;clc;
ITER=100; % number of random


addpath(genpath('~/toolboxes/'))
if ismac
    base_corpus='~/ResearchMIT/mixtures/timit-train/';
    addpath('~/ResearchMIT/toolboxes/Sound_Texture_Synthesis_Toolbox/');
else
    base_corpus='~/data/sounds2/timit-train/';
    addpath('~/mixtures');
    
    base_corpus='~/data/sounds/audiobooks-alex-wavs'
    
    
end
base_ptrn='s*.wav';

feature_fname='~/data/mixture-res/FEATURES-timit-mask-mix-SPEC-WHATTYPE-1-100.mat';

MYRMS=0; % rms of the mixture!!! important
%MYRMS=80; % rms of the mixture!!! important
is_sound=false;
is_show_snap=false; %for debuging showing snapshots
is_save_sound_with_features=true; %saving the real audio

%WHATTYPE=1; %normal mix
WHATTYPE=2; %ideal mask mix vs previous mix
%WHATTYPE=3; %ideal mask mix vs smooth gaussian mix

if WHATTYPE==1  
    IS_STAY_MASK=false;
    IS_NORM_MIX=true;
elseif WHATTYPE==2
    IS_STAY_MASK=true;
    IS_NORM_MIX=false;
else
    IS_STAY_MASK=false;
    IS_NORM_MIX=false;
end

%is_show_snap=true; %for debuging showing snapshots
%IS_STAY_MASK=true;
%IS_NORM_MIX=true;

MN=2;
FEATURES=cell(MN,1);
INFO=cell(MN,1);

FS0=16000;
DUR=1; %sec
M=2;

%NFEAT=1280;  % for 40*32 100 msec env subbands
NFEAT=640;  % for modpower
% NFEAT=1024;   %for C1 and C
% NFEAT=32*6*2; % for C2 (384)
% NFEAT=8288; %for all stat features
% NFEAT=32;  % for env_mean
%%

fprintf('reading corpus...\n');

cd (base_corpus);
files=dir(base_ptrn);
NF=length(files);
tic
fprintf('finding low volume threshold...\n');

VOLUMEITER=100;percent_remove=.18;
MYTHRESH=compute_rms(files,VOLUMEITER,percent_remove,DUR);
toc;


fprintf('making corpus...\n');
tic
mlength=inf;

for mm=1:2,
    
    
    
    myFEATURES=nan(ITER,NFEAT);
    myINFO=cell(ITER,1);
    
    
    cd (base_corpus);
    
    for KK=1:ITER
        
        fprintf('now in iteration %d of %d\t block number M=%d\n',KK,ITER,mm);
        mixme=cell(M,1);
        mynames=cell(M,1);
        
        cd (base_corpus);
        
        for I=1:M,
            
            myrms=0;
            while myrms<MYTHRESH
                
                smpls=0;fs=0;
                while (smpls-fs*DUR)<=0
                    iD=randi(NF,1,1);
                    
                    fname=files(iD).name;
                    info=audioinfo(fname);
                    smpls=info.TotalSamples;
                    fs=info.SampleRate;
                    if (smpls-fs*DUR)<0
                        fprintf('...too short \n')
                    end
                    
                end
                mypos=randi(smpls-fs*DUR,1,1);
                myrange=[mypos, mypos+fs*DUR];
                if myrms>0
                    fprintf('had to skip %d because rms was too small\n',KK);
                end
                
                [Y, FS]=audioread(fname, myrange);
                if size(Y,2)==2
                    Y=sum(Y,2);
                end
                myrms=sqrt(mean(Y.^2));
                
            end
            
            Y=double(Y);
            Y=0.03*(Y/sqrt(mean((Y.^2))));
            
            mixme{I}=resample(double(Y),FS0,FS);
            mixme{I}=0.03*mixme{I}/sqrt(mean(mixme{I}.^2));
            
            mixme{I}=mixme{I}*(10^(-MYRMS*(I-1)/20)); %note we make sounds uneven!
            
            mynames{I}=fname;
            
        end
        mlength=DUR;
        
        for I=1:M,
            mixme{I}=mixme{I}(1:round(FS0*mlength));
        end
        soundM=zeros(round(FS0*mlength),1);
        for I=1:M,
            soundM(1:round(FS0*mlength))=soundM(1:round(FS0*mlength))+ mixme{I}(1:round(FS0*mlength));
        end
        soundM=0.03*(soundM/sqrt(mean((soundM.^2))));
        %maskM=sqrt((temp_1))>(10^(MYRMS/20))*sqrt((temp_2));
        
        if IS_NORM_MIX
            if mm==1
                myts=mixme{1};
                myfs=FS0;
            else
                myts=soundM;
                myfs=FS0;
            end
            
        else
  
            S0 = nori_generate_fqsubands(soundM, FS0);
            if (mm==1)
                S1 = nori_generate_fqsubands(mixme{1}, FS0);
                S2 = nori_generate_fqsubands(mixme{2}, FS0);
                maskM=S1.subband_envs>S2.subband_envs;
                %maskM=20*log10(S1.subband_envs)>-15;
            end
            
            if (mm==2) && IS_STAY_MASK
                S1 = nori_generate_fqsubands(mixme{1}, FS0);
                S2 = nori_generate_fqsubands(mixme{2}, FS0);
                maskM=S1.subband_envs>S2.subband_envs;
                
                if KK>1
                    mask_old2=maskM;
                    maskM=mask_old;
                    mask_old=mask_old2;
                else
                    mask_old=maskM;
                    maskM=conv2(randn(size(S0.subband_envs)),ones(20),'same')>0;
                end
            end
            
            
            if mm==2 && (~IS_STAY_MASK)
                maskM=conv2(randn(size(S0.subband_envs)),ones(20),'same')>0;
            end
            
            
            
            collapse_subband_envs_1=S0.subband_envs.*maskM;
            collapse_subband_envs_2=S0.subband_envs.*(~maskM);
            
            collapse_subband_envs_1_up=resample(collapse_subband_envs_1,S0.P.audio_sr,S0.P.env_sr);
            collapse_subbands_1=S0.subbands.*collapse_subband_envs_1_up;
            collapse_sound_1=collapse_subbands(collapse_subbands_1,S0.audio_filts);
            
            collapse_subband_envs_2_up=resample(collapse_subband_envs_2,S0.P.audio_sr,S0.P.env_sr);
            collapse_subbands_2=S0.subbands.*collapse_subband_envs_2_up;
            collapse_sound_2=collapse_subbands(collapse_subbands_2,S0.audio_filts);
            
            sound_0=soundM;
            %else
            
            
            
            
            
            
            mfname=sprintf('mix-v1.M.%d.%d.wav',mm,KK);
            
           
            
            myts=collapse_sound_1;
            myfs=FS0;
        end
        
        
        
         
            
        tt3=1:length(myts);
        myts=myts/sqrt(mean(myts.^2))*0.03; %norma lize wav with rms
        
        
        S = nori_measure_texture_stats(myts, myfs); %measure stats
        S=S.S;
        
        feature=reshape(S.mod_power,[1 (length(S.Hz_mod_cfreqs)*length(S.audio_cutoffs_Hz))]); % rehshape stats
        myINFO{KK,1}.xlgnd=S.Hz_mod_cfreqs;
        myINFO{KK,1}.ylgnd=S.audio_cutoffs_Hz;
        myINFO{KK,1}.xlgnd_name='mod-fq(Hz)';
        myINFO{KK,1}.ylgnd_name='fq(Hz)';
        
%         
%         S = nori_generate_modsubands(myts, myfs);
%         
%         mymat=S.subband_envs';
%         feature=reshape(mymat,[1 numel(S.subband_envs)]); % rehshape stats
%         myINFO{KK,1}.xlgnd=1000*(1:size(mymat,2))/myfs;
%         myINFO{KK,1}.ylgnd=S.audio_cutoffs_Hz;
%         myINFO{KK,1}.xlgnd_name='tims(msec)';
%         myINFO{KK,1}.ylgnd_name='fq(Hz)';
%         
        
        
        
        
        
        %         if isempty(myFEATURES)
        %             NFEAT=length(feature);
        %             myFEATURES=nan(ITER,NFEAT);
        %         end
        %         assert(length(feature)==NFEAT)
        
        myFEATURES(KK,:)=feature;
        
        
        %myINFO{KK,1}.fname=mfname;
        %myINFO{KK,1}.range=[min(tt3),max(tt3)];
        %myINFO{KK,1}.Hz_mod_cfreqs=S.Hz_mod_cfreqs;
        %myINFO{KK,1}.audio_cutoffs_Hz=S.audio_cutoffs_Hz;
        
        if is_save_sound_with_features
         myINFO{KK,1}.audio=myts;
         myINFO{KK,1}.fs=myfs;
        end
        
        %downsampling
        %myINFO{KK,1}.audio=resample(myts,round(myfs/10),myfs);
        %myINFO{KK,1}.fs=round(myfs/10);
        
        if is_sound
                p1 = audioplayer(myts/max(myts), FS);p1.play
                pause(1.1*length(myts)/FS +0.5);
                
               
                figure(10);clf;
                %imagesc(S.Hz_mod_cfreqs,S.audio_cutoffs_Hz, S.mod_power);
                %clim=[min(min(S.mod_power)),max(max(S.mod_power))];
                %h=nori_log_imagesc(S.Hz_mod_cfreqs,S.audio_cutoffs_Hz, S.mod_power,clim,[]);axis xy;
                imagesc(1:32,(1:size(S.subband_envs,1))/S.env_sr, S.subband_envs');axis xy;
                
                drawnow
                %pause
                % audiowrite('MAR15_example1.wav',myts/max(myts),FS);
         end
        
    end
    FEATURES{mm}=myFEATURES;
    INFO{mm}=myINFO;
end

toc
%%


tic
fprintf('start saving to %s\n',feature_fname)
save(feature_fname,'-v7.3');
fprintf('end saving \n');
toc
%%
%if is_show_snap
figure(1);clf;imagesc([FEATURES{2};FEATURES{1}]);axis xy
%end
fprintf('DONE!');

%
%     %         jstats=struct();
%         %         [jstats.data,jstats.xleg,jstats.yleg,jstats.xlabels,jstats.ylabels,jstats.titles]=nori_stats_tocellarrays(S);
%         %         [jstats.vec,jstats.vecformat]=nori_cell_array_vectorize(jstats.data);
%         %
%         %         feature= jstats.vec;
%         %         myINFO{KK,1}.jstats=jstats;
%         %
%
%         feature=reshape(S.mod_power,[1 (length(S.Hz_mod_cfreqs)*length(S.audio_cutoffs_Hz))]); % rehshape stats
%         myINFO{KK,1}.xlgnd=S.Hz_mod_cfreqs;
%         myINFO{KK,1}.ylgnd=S.audio_cutoffs_Hz;
%         myINFO{KK,1}.xlgnd_name='mod-fq(Hz)';
%         myINFO{KK,1}.ylgnd_name='fq(Hz)';
%         %
%         %         feature=S.env_mean;
%         %         myINFO{KK,1}.xlgnd=S.audio_cutoffs_Hz;
%         %         myINFO{KK,1}.ylgnd=[1];
%         %         myINFO{KK,1}.xlgnd_name='fq(Hz)';
%         %         myINFO{KK,1}.ylgnd_name='1';
%         %
%         %         WC1=6;
%         %         feature=reshape(S.mod_C1(:,:,WC1),[1 (length(S.audio_cutoffs_Hz)*length(S.audio_cutoffs_Hz))]); % rehshape stats
%
%         %          feature=reshape(S.env_C(:,:),[1 (length(S.audio_cutoffs_Hz)*length(S.audio_cutoffs_Hz))]); % rehshape stats
%         %         myINFO{KK,1}.xlgnd=S.audio_cutoffs_Hz;
%         %         myINFO{KK,1}.ylgnd=S.audio_cutoffs_Hz;
%         %         myINFO{KK,1}.xlgnd_name='Coch. Channel (Hz)[C]';
%         %         myINFO{KK,1}.ylgnd_name='Coch. Channel (Hz)[C]';
%
%         %feature=reshape([S.mod_C2(:,:,1),S.mod_C2(:,:,2)],[1 32*6*2]); % rehshape stats
%         %         myINFO{KK,1}.xlgnd=1:12;
%         %         myINFO{KK,1}.ylgnd=S.audio_cutoffs_Hz;
%         %         myINFO{KK,1}.xlgnd_name='Modulation Channel (#real),(#imag) ';
%         %         myINFO{KK,1}.ylgnd_name='Coch. Channel (Hz)[C]';
%