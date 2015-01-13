% compute features
close all;clc;clear all;
NFB=20;
Ms=[1:9];
MN=length(Ms);
ITER=1000;
FEATURES=cell(MN,1);
% LF=300;HF=3000;
LF=300;HF=2500;

for m=1:MN,
    M=Ms(m);
    myhist_all=zeros(NFB,NFB);
    for KK=1:ITER,
        mfname=sprintf('data/mix-v1.M.%d.%d.wav',M,KK);
        display(sprintf('%s',mfname));
        
        [ts,fs]=audioread(mfname);
        [FilterBank,SubbandFrequencies]=make_erb_cos_filters(length(ts),fs,NFB-2,LF,HF);
        Cgrm= generate_subbands(ts,FilterBank);
        myCgrm=abs(Cgrm);
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
        FEATURES{m}=[FEATURES{m};myhist(:)'];
        
        myhist_all=myhist_all+myhist;
        
    end
    figure(2);
    subplot(3,3,m);
    imagesc(SubbandFrequencies,SubbandFrequencies,myhist_all);
    title(sprintf('num spks %d',M));
end
save('FEATURES-flat-cxxxv1.mat');
