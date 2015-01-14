clear all;close all;clc;
ITER=1000; % number of random
Ms=[1:9]; %all possible number of speakers
addpath('~/mixtures');
base_corpus='~/mixtures/inst';
output_data='~/mixtures/datainst';
  

length_min=1.5;
FS0=16000;
BEG=2;

    disp('reading corpus...');
    cd (base_corpus);
    files=dir('*.aiff');

NF=length(files);

for mm=1:length(Ms),
    M=Ms(mm);
    
    cd (output_data);
    for KK=1:ITER,
        display(sprintf('creating file %d of %d\t M=%d',KK,ITER,M));
        mixme=cell(M,1);
        mynames=cell(M,1);
        cd (base_corpus);
        mlength=length_min;
        for I=1:M,
            iD=randi(NF,1,1);
            
            fname=files(iD).name
%         [Y,FS,NBITS,CHUNKDATA] = aiffread(fname,[BEG*44100,(BEG+1)*44100]);
        [Y,FS,NBITS,CHUNKDATA] = aiffread(fname);
        assert(FS==44100);
        Y=double(Y);
        Y=Y/max(Y);
        Y=0.03*(Y/sqrt(mean((Y.^2))));
        R=resample(Y,FS0,FS);
        
        [~,idx]=max(abs(R));
        idx=min(max(idx,1+length_min*FS0/2),length(R)-length_min*FS0/2 -1);
        R=R(round((idx-length_min*FS0/2):(idx+length_min*FS0/2)));
        lengthR=size(R,1)/FS0;
            
            
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
        cd (output_data);
            wavwrite(soundM,FS0,mfname);
    end
end

% AP=audioplayer(soundM,FS);
% AP.play






