% compute features
%cd 'C:\Users\user\Dropbox (PPCA)\Research MIT\mixtures'
close all;clc;clear all;
NFB=150;
Ms=[2];
MN=length(Ms);
ITER=1000;
FEATURES=cell(MN,1);
% LF=300;HF=3000;
% LF=300;HF=2500;
LF=40;HF=8000;

figure(1);clf;
for m=1:MN,
    M=Ms(m);
    myhist_all=zeros(NFB,NFB);
    for KK=1:ITER,
        mfname=sprintf('data/mix-v1.M.%d.%d.wav',M,KK);
        display(sprintf('%s',mfname));
        
        [ts,fs]=audioread(mfname);
        ts=ts/sqrt(mean(ts.*ts));
        [FilterBank,SubbandFrequencies]=make_erb_cos_filters(length(ts),fs,NFB-2,LF,HF);
        Cgrm= generate_subbands(ts,FilterBank);
        myCgrm=sqrt((Cgrm.*Cgrm));
        %         myCgrm=(Cgrm);
        myhist=zeros(NFB,NFB);
        JUMP=50;
        for I=1:NFB,
            for J=1:NFB,
                %                  if I~=J
                myhist(I,J)=sum( (myCgrm(1:(end-JUMP),I)).*( myCgrm((1+JUMP):end,J)))/length(myCgrm(:,J));
                %                  else
                %                      myhist(I,J)=sum( (myCgrm(:,I)).*( myCgrm(:,J)))/length(myCgrm(:,J));
                %                      myhist(I,J)=sum( (myCgrm(:,I)).*( myCgrm(:,J)))/(10*length(myCgrm(:,J)));
                %                  end
            end
        end
        
        myhist=myhist/sum(sum(myhist));
        marg=sum(myhist);
        myhist_marg=repmat(marg,size(myhist,1),1).*repmat(marg',1,size(myhist,1));
        pmi=log2(myhist+eps)./(log2(myhist_marg+eps));
        
        subplot(2,3,1);
        imagesc(SubbandFrequencies,SubbandFrequencies,log2(myhist+eps));title('hist');axis xy;
        colorbar;
        subplot(2,3,2);
        imagesc(SubbandFrequencies,SubbandFrequencies,log2(myhist_marg+eps));title('hist-marg');axis xy;
        colorbar;
        
        subplot(2,3,3);
        imagesc(SubbandFrequencies,SubbandFrequencies,pmi);title('pmi');axis xy;
        colorbar;
        subplot(2,3,4);
        imagesc(1000*(1:size(myCgrm,1))/fs,SubbandFrequencies,(myCgrm'+eps));title('spectrogram');axis xy;
        
        subplot(2,3,5);
        imagesc(1000*(1:size(myCgrm,1))/fs,SubbandFrequencies,log2(myCgrm'+eps));title('spectrogram');axis xy;
        
        
        
        FEATURES{m}=[FEATURES{m};myhist(:)'];
        %         imagesc(SubbandFrequencies,SubbandFrequencies,myhist);
        figure(1);
        
        
        myhist_all=myhist_all+myhist;
        [C, ~, ~] = SpectralClustering(pmi, 2, 2);
        C=full(C);
        
        pos1=C(:,1)==1;
        pos2=C(:,2)==1;
        
        
        Cgrm1=zeros(size(Cgrm));Cgrm1(:,pos1)=Cgrm(:,pos1);
        Cgrm2=zeros(size(Cgrm));Cgrm2(:,pos2)=Cgrm(:,pos2);
        
        
        %If you want to generate a time series from a Cochleagram use
        Y1=0.5*collapse_subbands(Cgrm1,FilterBank);
        Y2=0.5*collapse_subbands(Cgrm2,FilterBank);
        drawnow;
        
        for tt=1:2,
        
        p = audioplayer(0.1*ts, fs);p.play
        pause(1);
        p1 = audioplayer(0.1*Y1, fs);p1.play
        pause(1);
        p2 = audioplayer(0.1*Y2, fs);p2.play
        pause(1);
        end
        
        
        subplot(2,3,6);
        clust=zeros(size(myCgrm'));
        clust(C(:,1)==1,:)=1;
        clust(C(:,2)==1,:)=2;
        imagesc(1000*(1:size(myCgrm,1))/fs,SubbandFrequencies,clust);title('clustered-spectrogram');axis xy;
        %         imagesc(1000*(1:size(myCgrm,1))/fs,SubbandFrequencies,(myCgrm'+eps));title('spectrogram');axis xy;
        
        
        
        pause;
        
    end
    figure(2);
    subplot(3,3,m);
    imagesc(SubbandFrequencies,SubbandFrequencies,myhist_all);
    title(sprintf('num spks %d',M));
end
save('FEATURES-flat-cxxxv1.mat');
