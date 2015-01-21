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
%START=70
START=3;
for m=1:MN,
    M=Ms(m);
    myhist_all=zeros(NFB,NFB);
    for KK=START:ITER,
        mfname=sprintf('data-inst/mix-v1.M.%d.%d.wav',M,KK);
        display(sprintf('%s',mfname));
        
        %[ts,fs]=audioread(mfname,[0.2*16000 1.4*16000]);
        [ts,fs]=audioread(mfname);
        
        ts=ts/sqrt(mean(ts.*ts));
        [FilterBank,SubbandFrequencies]=make_erb_cos_filters(length(ts),fs,NFB-2,LF,HF);
        Cgrm= generate_subbands(ts,FilterBank);
        myCgrm=sqrt((Cgrm.*Cgrm));
        %         myCgrm=(Cgrm);
        myhist=zeros(NFB,NFB);
        JUMP= 0;
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
%         subplot(2,3,2);
%         imagesc(SubbandFrequencies,SubbandFrequencies,log2(myhist_marg+eps));title('hist-marg');axis xy;
%         colorbar;
    
        
        subplot(2,3,3);
        %imagesc(SubbandFrequencies,SubbandFrequencies,pmi);title('pmi');axis xy;
        %imagesc(myhist)
        nm=(myhist+eps);
        %nm=nm./repmat(sum(nm),[size(nm,1),1]);
        
        imagesc(nm);title('nm');axis xy;
        
        colorbar;
        subplot(2,3,4);
         imagesc(1000*(1:size(myCgrm,1))/fs,SubbandFrequencies,log2(abs((myCgrm'+eps))));title('spectrogram');axis xy;
         
     
        FEATURES{m}=[FEATURES{m};myhist(:)'];
        %         imagesc(SubbandFrequencies,SubbandFrequencies,myhist);
        figure(1);
        
        
        myhist_all=myhist_all+myhist;
  %       [C, ~, ~] = SpectralClustering(log2(pmi+eps), 2, 2);
    %    [C, ~, ~] = SpectralClustering(log2(myhist), 2, 2);
           [C, ~, ~] = SpectralClustering(nm, 2, 2);
%         sigma = .1;
% nbclusters = 2;
% 
% [C, ~,~] = spcl(log2(myhist)/max(log2(myhist)), nbclusters, sigma,'np' , [2 2]);
% 
% C=[C==1,C==2];  
%         [pc, net] = nlpca(myCgrm(1:10:end), k)
% 
% pc = nlpca_get_components(net, data)
% data_reconstruction = nlpca_get_data(net, pc)
      
        
        
        C=full(C);
        
        pos1=C(:,1)==1;
        pos2=C(:,2)==1;
        pos3=[find(pos1);find(pos2)];
        
        Cgrm1=zeros(size(Cgrm));Cgrm1(:,pos1)=Cgrm(:,pos1);
        Cgrm2=zeros(size(Cgrm));Cgrm2(:,pos2)=Cgrm(:,pos2);
        
        subplot(2,3,2);
        imagesc(log2(myhist(pos3,pos3)+eps));title('re-ordered hist');axis xy;hold on;
        plot([sum(pos1),sum(pos1)],[0,length(pos3)],'k-')
        plot([0,length(pos3)],[sum(pos1),sum(pos1)],'k-')
        colorbar;
        
        
        %If you want to generate a time series from a Cochleagram use
        Y1=0.5*collapse_subbands(Cgrm1,FilterBank);
        Y2=0.5*collapse_subbands(Cgrm2,FilterBank);
        
        %%
%        [pc, net] = nlpca(log2(myCgrm(1:10:end,:)+eps)', 2)
%        %%
%        pc = nlpca_get_components(net, data)
% % data_reconstruction = nlpca_get_data(net, pc)
%         %%
%         data_reconstruction1 = nlpca_get_data(net, [pc(1,:);0*pc(2,:)]);
%         data_reconstruction2 = nlpca_get_data(net, [0*pc(1,:);pc(2,:)]);
%         %data_reconstruction1=resample(data_reconstruction1,10,1);
%         %data_reconstruction2=resample(data_reconstruction2,10,1);
%         %data_reconstruction1=data_reconstruction1(size(myCgrm(1:10:end,:)));
%         %data_reconstruction2=data_reconstruction2(size(myCgrm(1:10:end,:)));
%         
%         %%
%         figure(100);clf;imagesc(log2(myCgrm(1:10:end,:)'+eps))
%         figure(101);clf;imagesc(log2(abs(data_reconstruction1+eps)))
%         figure(102);clf;imagesc(log2(abs(data_reconstruction2+eps)))
%         %%
%         drawnow;
%        p = audioplayer(0.1*ts, fs);p.play
%         pause(1.5);
%         my1=Cgrm./sqrt(abs(Cgrm.*Cgrm)).*(2.^data_reconstruction1);
%         mY1=0.5*collapse_subbands(my1,FilterBank);
%         p = audioplayer(mY1, fs);p.play
%         pause(1);
%         
%         my2=Cgrm./sqrt(abs(Cgrm.*Cgrm)).*(2.^data_reconstruction2);
%         mY2=0.5*collapse_subbands(my2,FilterBank);
%         p = audioplayer(mY2, fs);p.play
        
        %%
        Cgrm1p=Cgrm1+0.0051*Cgrm2;
        Cgrm2p=Cgrm2+0.0051*Cgrm1;
           subplot(2,3,5);
        imagesc(1000*(1:size(myCgrm,1))/fs,SubbandFrequencies,log2(sqrt((Cgrm1p.*Cgrm1p))'+eps));title('spectrogram1');axis xy;
        
        subplot(2,3,6);
        imagesc(1000*(1:size(myCgrm,1))/fs,SubbandFrequencies,log2(sqrt((Cgrm2p.*Cgrm2p))'+eps));title('spectrogram2');axis xy;
        
        
        drawnow;
        
        for tt=1:2,
        
        p = audioplayer(0.1*ts, fs);p.play
        pause(1.5*length(ts)/fs);
        p1 = audioplayer(0.1*Y1, fs);p1.play
        pause(1.1*length(ts)/fs);
        p2 = audioplayer(0.1*Y2, fs);p2.play
        pause(1.1*length(ts)/fs);
        pause
        end
        
        
%         subplot(2,3,6);
%         clust=zeros(size(myCgrm'));
%         clust(C(:,1)==1,:)=1;
%         clust(C(:,2)==1,:)=2;
%         imagesc(1000*(1:size(myCgrm,1))/fs,SubbandFrequencies,clust);title('clustered-spectrogram');axis xy;
        %         imagesc(1000*(1:size(myCgrm,1))/fs,SubbandFrequencies,(myCgrm'+eps));title('spectrogram');axis xy;
        
        
        
        pause;
        
    end
    figure(2);
    subplot(3,3,m);
    imagesc(SubbandFrequencies,SubbandFrequencies,myhist_all);
    title(sprintf('num spks %d',M));
end
save('FEATURES-flat-cxxxv1.mat');
