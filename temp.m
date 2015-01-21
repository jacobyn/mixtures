clear all;close all;clc;
ITER=10; % number of random
Ms=[2]; %all possible number of speakers
%addpath('C:\Users\user\Dropbox (PPCA)\Research MIT\mixtures');
addpath('~/ResearchMIT/mixtures/');
addpath(genpath('~/ResearchMIT/toolboxes'));
base_corpus='~/ResearchMIT/mixtures/timit-train';
%output_data='C:\Users\user\Dropbox (PPCA)\Research MIT\mixtures\data-solopiece';
LF=40;HF=8000;NFB=150;

FS0=16000;

disp('reading corpus...');
cd (base_corpus);
files=dir('s*.wav');
%%
NF=length(files);

for mm=1:length(Ms),
    M=Ms(mm);
    
    
    for KK=1:ITER,
        
        display(sprintf('now in file %d of %d\t M=%d',KK,ITER,M));
        mixme=cell(M,1);
        myCgrm=cell(M,1);
        mynames=cell(M,1);
        cd (base_corpus);
        mlength=inf;
        for I=1:M,
            iD=randi(NF,1,1);
            
            fname=files(iD).name;
            %             info=audioinfo(fname);
            %             smpls=info.TotalSamples;
            %             fs=info.SampleRate;
            %             assert(smpls-fs*DUR>0);
            %             mypos=randi(smpls-fs*DUR,1,1);
            %             myrange=[mypos, mypos+fs*DUR];
            
            [Y, FS]=audioread(fname);
            
            %         [Y,FS,NBITS,CHUNKDATA] = aiffread(fname,[BEG*44100,(BEG+1)*44100]);%             [Y,FS,NBITS,CHUNKDATA] = aiffread(fname);
            
            %assert(FS==44100);
            Y=double(Y);
            Y=0.03*(Y/sqrt(mean((Y.^2))));
            R=resample(Y,FS0,FS);
            lengthR=length(R)/FS0;
            mixme{I}=0.03*(R/sqrt(mean((R.^2))));
            mlength=min(mlength,lengthR);
            mynames{I}=fname;
        end
        
        soundM=zeros(round(FS0*mlength),1);
        for I=1:M,
            mixme{I}=mixme{I}(1:round(FS0*mlength));
            soundM=soundM+mixme{I};
        end
        soundM=0.03*(soundM/sqrt(mean((soundM.^2))));
        
        [FilterBank,SubbandFrequencies]=make_erb_cos_filters(length(soundM),FS,NFB-2,LF,HF);
        myCgrmM= generate_subbands(soundM,FilterBank);
        for I=1:M,
            [FilterBank,SubbandFrequencies]=make_erb_cos_filters(length(mixme{I}),FS,NFB-2,LF,HF);
            myCgrm{I}= generate_subbands(mixme{I},FilterBank);
        end
        
        assert(M==2);
        maskM=(myCgrm{1}.^2)>(myCgrm{2}.^2);
        Y1=collapse_subbands(myCgrmM.*(maskM==0),FilterBank);
        Y2=collapse_subbands(myCgrmM.*(maskM==1),FilterBank);
    
        
        figure(1);
        subplot(3,2,1);imagesc(log(myCgrmM.^2)');axis xy;title('mixture');
        subplot(3,2,2);imagesc(maskM');axis xy; title('mask');
        subplot(3,2,3);imagesc(log(myCgrm{1}.^2)');axis xy; title('original 1');
        subplot(3,2,4);imagesc(log(myCgrm{2}.^2)');axis xy; title('original 2');
        
        subplot(3,2,5);imagesc(log((myCgrm{1}.*(maskM==0)).^2 )');axis xy;title('masked mixture 1');
        subplot(3,2,6);imagesc(log((myCgrm{2}.*(maskM==1)).^2 )');axis xy;title('masked mixture 2');
       
        p = audioplayer(soundM, FS);p.play
        pause(1.5*length(soundM)/FS);
        %          p = audioplayer(mixme{1}, FS);p.play
        %          pause(1.1*length(mixme{1})/FS);
        %          p = audioplayer(mixme{2}, FS);p.play
        %          pause(1.1*length(mixme{2})/FS);
        
        
        
        p1 = audioplayer(Y1, FS);p1.play
        pause(1.1*length(Y1)/FS);
        p2 = audioplayer(Y2, FS);p2.play
        pause(1.1*length(Y2)/FS);
        
        
        %%
        %         mfname=sprintf('mix-v1.M.%d.%d.wav',M,KK);
        
        %
        %         cd (output_data);
        %         %wavwrite(soundM,FS0,mfname);
        %         audiowrite(mfname,soundM,FS0)
        
    end
end







