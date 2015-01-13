% compute features
close all;clc;clear all;
is_flat=false;
NFB=50;
Ms=[1,9];
MN=length(Ms);
ITER=1000;
FEATURES=cell(MN,1);
% LF=300;HF=3000;
LF=300;HF=2500;
% LF=40;HF=10000;

for m=1:MN,
    M=Ms(m);
    if ~is_flat
        myhist_all=zeros(NFB,NFB);
    else
        myhist_all=zeros(1,NFB);
    end
    for KK=1:ITER,
        mfname=sprintf('data/mix-v1.M.%d.%d.wav',M,KK);
        display(sprintf('%s',mfname));
        
        [ts,fs]=audioread(mfname);
        ts=ts/max(abs(ts));
        [FilterBank,SubbandFrequencies]=make_erb_cos_filters(length(ts),fs,NFB-2,LF,HF);
        Cgrm= generate_subbands(ts,FilterBank);
        myCgrm=abs(Cgrm);
        imagesc(1000*(1:size(myCgrm,1))/fs,SubbandFrequencies,myCgrm');axis xy;
        
        pause
        
        if ~is_flat
            myhist=zeros(NFB,NFB);
            for I=1:NFB,
                for J=1:NFB,
                    if I~=J
                        myhist(I,J)=sum( (myCgrm(:,I)).*( myCgrm(:,J)))/length(myCgrm(:,J));
                    else
                        myhist(I,J)=sum( (myCgrm(:,I)).*( myCgrm(:,J)))/(10*length(myCgrm(:,J)));
                    end
                end
            end
            myhist=myhist/max(max(myhist));
        else
            myhist=sum(myCgrm.*myCgrm,1);
            myhist=myhist/max(max(myhist));
        end
        
        
        FEATURES{m}=[FEATURES{m};myhist(:)'];
        
        myhist_all=myhist_all+myhist;
        
    end
    figure(2);
    subplot(3,3,m);
    imagesc(SubbandFrequencies,SubbandFrequencies,myhist_all);
    title(sprintf('num spks %d',M));
end
save('FEATURES-v7x.mat');
