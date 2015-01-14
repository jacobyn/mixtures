clear all;close all;clc;
ITER=1000; % number of random
Ms=[1:2]; %all possible number of speakers
addpath('C:\Users\user\Dropbox (PPCA)\Research MIT\mixtures');
base_corpus='C:\Users\user\Dropbox (PPCA)\Research MIT\solo-pieces';
output_data='C:\Users\user\Dropbox (PPCA)\Research MIT\mixtures\data-solopiece';


FS0=16000;
DUR=1.0; %sec

disp('reading corpus...');
cd (base_corpus);
files=dir('*.mp3');
%%
NF=length(files);
mlength=inf;
for mm=1:length(Ms),
    M=Ms(mm);
    
    cd (output_data);
    for KK=1:ITER,
%         if (M<2) && (KK<350)
%             fprintf('skipping....\n')
%             continue %DELMED!!!!!
%         end
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
            
            %         [Y,FS,NBITS,CHUNKDATA] = aiffread(fname,[BEG*44100,(BEG+1)*44100]);
%             [Y,FS,NBITS,CHUNKDATA] = aiffread(fname);
            assert(FS==44100);
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
        
        
        cd (output_data);
        %wavwrite(soundM,FS0,mfname);
        audiowrite(mfname,soundM,FS0)
        
    end
end







