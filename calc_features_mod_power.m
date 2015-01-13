
% compute features
close all;clc;clear all;
Ms=[1,2]; %num of spks in mixture Ms are vector of all m=#spk to compute
ITER=1000;

% datadir='data-inst\'; %instruments
datadir='data\'; %speech
datadir='C:\Users\user\Dropbox (PPCA)\Research MIT\mixtures\data-solopiece\'; %solo cello



outdir='C:\Users\user\Dropbox (PPCA)\Research MIT\mixtures\';
cd (outdir);
addpath('C:\Users\user\Dropbox (PPCA)\Research MIT\Sound_Texture_Synthesis_Toolbox');
addpath('C:\Users\user\Dropbox (PPCA)\Research MIT\mixtures\');


is_show_snap=false; %for debuging showing snapshots
% is_show_snap=true; %for debuging showing snapshots

MN=length(Ms);
FEATURES=cell(MN,1);
INFO=cell(MN,1);


tsold=-1;
for m=1:MN,
    M=Ms(m);
    INFO{m}=cell(2,1);
    cnt=0;
    
    fnames=dir(sprintf('%s*.M.%d.*',datadir,m));
    for KK=1:ITER,
        mfname=sprintf('%s%s',datadir,fnames(KK).name);
        display(sprintf('%s',mfname));
        
        [ts,fs]=audioread(mfname);
        
        tt3=1:length(ts);
        ts=ts/sqrt(mean(ts.^2))*0.03; %normalize wav with rms
        
%         ll=floor(length(ts)/2)*2; % Josh requires a muliple of 2!! (correct that inside Josh's lib, this should not be like that)
%         S = nori_measure_texture_stats(ts(1:ll), fs); %measure stats        
        S = nori_measure_texture_stats(ts, fs); %measure stats        
        feature=reshape(S.mod_power,[1 (length(S.Hz_mod_cfreqs)*length(S.audio_cutoffs_Hz))]); % rehshape stats
        
        if is_show_snap
            imagesc(S.Hz_mod_cfreqs,S.audio_cutoffs_Hz, S.mod_power);axis xy;
            title(fnames(KK).name);  
            PP=audioplayer(ts,fs);PP.play;
            pause(1.1*length(ts)/fs);
        end
        cnt=cnt+1;
        FEATURES{m}=[FEATURES{m};feature];
        INFO{m}{cnt,1}.fname=mfname;
        INFO{m}{cnt,1}.range=[min(tt3),max(tt3)];
        INFO{m}{cnt,1}.Hz_mod_cfreqs=S.Hz_mod_cfreqs;
        INFO{m}{cnt,1}.audio_cutoffs_Hz=S.audio_cutoffs_Hz;
        
        
    end
end
cd (outdir);
% save('FEATURES-N-v4.mat');
save('FEATURES-cello-v2.mat');

figure(1);clf;imagesc([FEATURES{2};FEATURES{1}]);axis xy
% figure(1);imagesc([log(abs(FEATURES{2}));log(abs(FEATURES{1}))]);axis xy