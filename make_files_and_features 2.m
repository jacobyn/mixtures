%matlabpool
clear all;close all;clc;
ITER=1000; % number of random
Ms=[1:2]; %all possible number of speakers


base_corpus='~/mcdermott/shared/Audiobooks';
base_ptrn='Hun*.wav';
base_corpus='~/mixtures/data/music-cello';
base_ptrn='flu*wav';
base_ptrn='*wav';


output_audiodata='~/mixtures/data/audiobooksall';
% feature_fname='~/mixtures/FEATURES-AUDIO-BOOKS-HUN.mat';
feature_fname='~/mixtures/FEATURES-3INST.mat';

addpath('~/Sound_Texture_Synthesis_Toolbox');
addpath('~/mixtures/');


is_show_snap=false; %for debuging showing snapshots
% is_show_snap=true; %for debuging showing snapshots


FS0=16000;
DUR=1.0; %sec
VOLUMEITER=100;
fprintf('reading corpus...\n');
cd (base_corpus);
files=dir(base_ptrn);
files
%%
fprintf('making file\n');
tic
fprintf('computing threshold for rms computation\n');
NF=length(files);

rmss=nan(VOLUMEITER,1);

for KK=1:VOLUMEITER,
    tiD=randi(NF,1,1);
    
    tfname=files(tiD).name;
    tinfo=audioinfo(tfname);
    tsmpls=tinfo.TotalSamples;
    tfs=tinfo.SampleRate;
    assert(tsmpls-tfs*DUR>0);
    tmypos=randi(tsmpls-tfs*DUR,1,1);
    tmyrange=[tmypos, tmypos+tfs*DUR];
    
    [tY, tFS]=audioread(tfname, tmyrange);
    if size(tY,2)==2
        tY=sum(tY,2);
    end
    rmss(KK)=sqrt(mean(tY.^2));
    
end
srms=sort(rmss);
MYTHRESH=srms(round(VOLUMEITER*.18));
toc
tic
NF=length(files);

mlength=inf;
for mm=1:length(Ms),
    M=Ms(mm);
    
    cd (output_audiodata);
    
    
    cd (base_corpus);
    
    parfor (KK=1:ITER)
        
        display(sprintf('creating file %d of %d\t M=%d',KK,ITER,M));
        mixme=cell(M,1);
        mynames=cell(M,1);
        
        cd (base_corpus);

        for I=1:M,
            
            myrms=0;
            while myrms<MYTHRESH
                iD=randi(NF,1,1);
                fname=files(iD).name;
                info=audioinfo(fname);
                smpls=info.TotalSamples;
                fs=info.SampleRate;
                assert(smpls-fs*DUR>0);
                mypos=randi(smpls-fs*DUR,1,1);
                myrange=[mypos, mypos+fs*DUR];
                if myrms>0
                    fprintf('had to skip %d\n',KK);
                end
                
                [Y, FS]=audioread(fname, myrange);
                if size(Y,2)==2
                    Y=sum(Y,2);
                end
                myrms=sqrt(mean(Y.^2));
                
            end
            %%% note note! no normalization
            %             Y=double(Y);
            %             Y=Y/max(Y);
            %             Y=0.03*(Y/sqrt(mean((Y.^2))));
            %             R=resample(Y,FS0,FS);
            %             lengthR=length(R)/FS0;
            %             mixme{I}=0.03*(R/sqrt(mean((R.^2))));
            
            mixme{I}=resample(double(Y),FS0,FS);
            
            
            mynames{I}=fname;
        end
        mlength=DUR;
        soundM=zeros(round(FS0*mlength),1);
        for I=1:M,
            soundM(1:round(FS0*mlength))=soundM(1:round(FS0*mlength))+ mixme{I}(1:round(FS0*mlength));
        end
        %         soundM=0.03*(soundM/sqrt(mean((soundM.2))));
        soundM=0.03*soundM;
        mfname=sprintf('mix-v1.M.%d.%d.wav',M,KK);
        
        
        cd (output_audiodata);
        %wavwrite(soundM,FS0,mfname);
        audiowrite(mfname,soundM,FS0)
        
    end
end
toc
%%
fprintf('makinmg feature\n');
tic
MN=length(Ms);
FEATURES=cell(MN,1);
INFO=cell(MN,1);
NFEAT=640;

cd (output_audiodata);
tsold=-1;
for m=1:MN,
    M=Ms(m);
    INFO{m}=cell(2,1);
    
    myFEATURES=nan(ITER,NFEAT);
    myINFO=cell(ITER,1);
    cnt=0;
    clear Y;
    clear FS;
    fnames=dir(sprintf('*.M.%d.*',m));
    parfor (KK=1:ITER)
        mfname=sprintf('%s',fnames(KK).name);
        display(sprintf('%s',mfname));
        
        [myts,myfs]=audioread(mfname);
        
        tt3=1:length(myts);
        myts=myts/sqrt(mean(myts.^2))*0.03; %normalize wav with rms
        
        S = nori_measure_texture_stats(myts, myfs); %measure stats
        feature=reshape(S.mod_power,[1 (length(S.Hz_mod_cfreqs)*length(S.audio_cutoffs_Hz))]); % rehshape stats
        
        if is_show_snap
            imagesc(S.Hz_mod_cfreqs,S.audio_cutoffs_Hz, S.mod_power);axis xy;
            title(fnames(KK).name);
            PP=audioplayer(ts,fs);PP.play;
            pause(1.1*length(ts)/fs);
        end
        assert(length(feature)==NFEAT)
        
        myFEATURES(KK,:)=feature;
        
        
        myINFO{KK,1}.fname=mfname;
        myINFO{KK,1}.range=[min(tt3),max(tt3)];
        myINFO{KK,1}.Hz_mod_cfreqs=S.Hz_mod_cfreqs;
        myINFO{KK,1}.audio_cutoffs_Hz=S.audio_cutoffs_Hz;
        myINFO{KK,1}.audio=myts;
        myINFO{KK,1}.fs=myfs;
        
        
    end
    FEATURES{m}=myFEATURES;
    INFO{m}=myINFO;
end
toc
save(feature_fname);
if is_show_snap
    figure(1);clf;imagesc([FEATURES{2};FEATURES{1}]);axis xy
end






