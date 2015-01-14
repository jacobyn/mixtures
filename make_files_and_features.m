clear all;close all;clc;
ITER=10; % number of random
Ms=[1:2]; %all possible number of speakers


base_corpus='C:\Users\user\Dropbox (PPCA)\Research MIT\audio_books';
base_ptrn='*.wav';
output_audiodata='C:\Users\user\Dropbox (PPCA)\Research MIT\mixtures\data-books';
feature_fname='C:\Users\user\Dropbox (PPCA)\Research MIT\mixtures\FEATURES-cello-v2.mat';

addpath('C:\Users\user\Dropbox (PPCA)\Research MIT\Sound_Texture_Synthesis_Toolbox');
addpath('C:\Users\user\Dropbox (PPCA)\Research MIT\mixtures\');


is_show_snap=false; %for debuging showing snapshots
% is_show_snap=true; %for debuging showing snapshots


FS0=16000;
DUR=1.0; %sec

fprintf('reading corpus...\n');
cd (base_corpus);
files=dir(base_ptrn);
%%
NF=length(files);
mlength=inf;
for mm=1:length(Ms),
    M=Ms(mm);
    
    cd (output_audiodata);
    for KK=1:ITER,
        
        display(sprintf('creating file %d of %d\t M=%d',KK,ITER,M));
        mixme=cell(M,1);
        mynames=cell(M,1);
        cd (base_corpus);
        
        for I=1:M,
            iD=randi(NF,1,1);
            
            fname=files(iD).name;
            info=audioinfo(fname);
            smpls=info.TotalSamples;
            fs=info.SampleRate;
            assert(smpls-fs*DUR>0);
            mypos=randi(smpls-fs*DUR,1,1);
            myrange=[mypos, mypos+fs*DUR];
            
            [Y, FS]=audioread(fname, myrange);
            
            Y=double(Y);
            Y=Y/max(Y);
            Y=0.03*(Y/sqrt(mean((Y.^2))));
            R=resample(Y,FS0,FS);
            lengthR=length(R)/FS0;
            mixme{I}=0.03*(R/sqrt(mean((R.^2))));
            mlength=min(mlength,lengthR);
            mynames{I}=fname;
        end
        
        soundM=zeros(round(FS0*mlength),1);
        for I=1:M,
            soundM(1:round(FS0*mlength))=soundM(1:round(FS0*mlength))+ mixme{I}(1:round(FS0*mlength));
        end
        soundM=0.03*(soundM/sqrt(mean((soundM.^2))));
        mfname=sprintf('mix-v1.M.%d.%d.wav',M,KK);
        
        
        cd (output_audiodata);
        %wavwrite(soundM,FS0,mfname);
        audiowrite(mfname,soundM,FS0)
        
    end
end
%%

MN=length(Ms);
FEATURES=cell(MN,1);
INFO=cell(MN,1);

cd (output_audiodata);
tsold=-1;
for m=1:MN,
    M=Ms(m);
    INFO{m}=cell(2,1);
    cnt=0;
    
    fnames=dir(sprintf('*.M.%d.*',m));
%     for KK=1:ITER,
    for KK=1:ITER,
        mfname=sprintf('%s',fnames(KK).name);
        display(sprintf('%s',mfname));
        
        [ts,fs]=audioread(mfname);
        
        tt3=1:length(ts);
        ts=ts/sqrt(mean(ts.^2))*0.03; %normalize wav with rms
        
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
save(feature_fname);
if is_show_snap
    figure(1);clf;imagesc([FEATURES{2};FEATURES{1}]);axis xy
end






