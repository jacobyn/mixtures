%matlabpoolv
clear all;close all;clc;
ITER=10000; % number of random
Ms=[1:2]; %all possible number of speakers


addpath(genpath('~/toolboxes/'))
base_corpus='~/data/sounds/timit-train/';
base_ptrn='s*.wav';
output_audiodata='~/data/mixture-temp';
feature_fname='~/data/mixture-res/FEATURES-timit-OK-10k.mat'


% addpath(genpath('~/toolboxes/'))
% base_corpus='~/data/audio-insts/';
% base_ptrn='cello*wav';
% output_audiodata='~/data/mixture-temp';
% feature_fname='~/data/mixture-res/FEATURES-cello-10k.mat'




%     base_ptrn='*.wav';
%     output_audiodata='C:\Users\user\Dropbox (PPCA)\Research MIT\mixtures\data-books';
%     feature_fname='C:\Users\user\Dropbox (PPCA)\Research MIT\mixtures\FEATURES-TEMP.mat';

% for mac:   addpath(genpath('~/ResearchMIT/toolboxes/'))

% if isunix
%     base_corpus='~/mcdermott/shared/Audiobooks';
%     base_ptrn='En*.wav';
%     %     base_corpus='~/mixtures/data/music-cello';
%     %      base_ptrn='flu*wav';
%     %     base_ptrn='pia*wav';
%     output_audiodata='~/mixtures/data/audiobooksall';
%     % feature_fname='~/mixtures/FEATURES-AUDIO-BOOKS-HUN.mat';
%     feature_fname='~/mixtures/FEATURES-ALLJSTATS-ENGLISH.mat'
%
%     addpath('~/Sound_Texture_Synthesis_Toolbox');
%     addpath('~/mixtures/');
% else
%     base_corpus='C:\Users\user\Dropbox (PPCA)\Research MIT\audio_books';
%     base_ptrn='*.wav';
%     output_audiodata='C:\Users\user\Dropbox (PPCA)\Research MIT\mixtures\data-books';
%     feature_fname='C:\Users\user\Dropbox (PPCA)\Research MIT\mixtures\FEATURES-TEMP.mat';
%     addpath('C:\Users\user\Dropbox (PPCA)\Research MIT\Sound_Texture_Synthesis_Toolbox');
%     addpath('C:\Users\user\Dropbox (PPCA)\Research MIT\mixtures\');
%
% end

NFEAT=640;  % for modpower
NFEAT=1024;   %for C1 and C
NFEAT=32*6*2; % for C2 (384)
NFEAT=8288; %for all stat features

is_show_snap=false; %for debuging showing snapshots
%  is_show_snap=true; %for debuging showing snapshots


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
    tsmpls=0;tfs=0;
    while (tsmpls-tfs*DUR)<=0
        tiD=randi(NF,1,1);
        
        tfname=files(tiD).name;
        tinfo=audioinfo(tfname);
        tsmpls=tinfo.TotalSamples;
        tfs=tinfo.SampleRate;
        if (tsmpls-tfs*DUR)<0
            fprintf('filename %s too short (%g sec)\n', tfname,tsmpls/tfs)
        end
        
    end
    
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
                
                smpls=0;fs=0;
                while (smpls-fs*DUR)<=0
                    iD=randi(NF,1,1);
                    
                    fname=files(iD).name;
                    info=audioinfo(fname);
                    smpls=info.TotalSamples;
                    fs=info.SampleRate;
                    if (smpls-fs*DUR)<0
                        fprintf('filename %s too short (%g sec)\n', tfname,tsmpls/tfs)
                    end
                    
                end
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
fprintf('making feature\n');
tic
MN=length(Ms);
FEATURES=cell(MN,1);
INFO=cell(MN,1);


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
        jstats=struct();
        [jstats.data,jstats.xleg,jstats.yleg,jstats.xlabels,jstats.ylabels,jstats.titles]=nori_stats_tocellarrays(S);
        [jstats.vec,jstats.vecformat]=nori_cell_array_vectorize(jstats.data);
        
        
        
        %         feature=reshape(S.mod_power,[1 (length(S.Hz_mod_cfreqs)*length(S.audio_cutoffs_Hz))]); % rehshape stats
        %         myINFO{KK,1}.xlgnd=S.Hz_mod_cfreqs;
        %         myINFO{KK,1}.ylgnd=S.audio_cutoffs_Hz;
        %         myINFO{KK,1}.xlgnd_name='mod_fq(Hz)';
        %         myINFO{KK,1}.ylgnd_name='fq(Hz)';
        %
        %         WC1=6;
        %         feature=reshape(S.mod_C1(:,:,WC1),[1 (length(S.audio_cutoffs_Hz)*length(S.audio_cutoffs_Hz))]); % rehshape stats
        
        %          feature=reshape(S.env_C(:,:),[1 (length(S.audio_cutoffs_Hz)*length(S.audio_cutoffs_Hz))]); % rehshape stats
        %         myINFO{KK,1}.xlgnd=S.audio_cutoffs_Hz;
        %         myINFO{KK,1}.ylgnd=S.audio_cutoffs_Hz;
        %         myINFO{KK,1}.xlgnd_name='Coch. Channel (Hz)[C]';
        %         myINFO{KK,1}.ylgnd_name='Coch. Channel (Hz)[C]';
        
        %feature=reshape([S.mod_C2(:,:,1),S.mod_C2(:,:,2)],[1 32*6*2]); % rehshape stats
        %         myINFO{KK,1}.xlgnd=1:12;
        %         myINFO{KK,1}.ylgnd=S.audio_cutoffs_Hz;
        %         myINFO{KK,1}.xlgnd_name='Modulation Channel (#real),(#imag) ';
        %         myINFO{KK,1}.ylgnd_name='Coch. Channel (Hz)[C]';
        feature= jstats.vec;
        myINFO{KK,1}.jstats=jstats;
        
        
        
        if is_show_snap
            imagesc(S.Hz_mod_cfreqs,S.audio_cutoffs_Hz, S.mod_power);axis xy;
            title(fnames(KK).name);
            PP=audioplayer(myts,myfs);PP.play;
            pause(1.1*length(myts)/myfs);
            pause
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
tic
fprintf('start saving to %s\n',feature_fname)
save(feature_fname,'-v7.3');
fprintf('end saving \n');
toc
if is_show_snap
    figure(1);clf;imagesc([FEATURES{2};FEATURES{1}]);axis xy
end
fprintf('DONE!');