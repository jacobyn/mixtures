clear all;close all;clc;
ITER=4000; % number of random
Ms=[2]; %all possible number of speakers
%addpath('C:\Users\user\Dropbox (PPCA)\Research MIT\mixtures');
addpath('~/ResearchMIT/mixtures/');
addpath(genpath('~/ResearchMIT/toolboxes'));
base_corpus='~/ResearchMIT/mixtures/timit-train';
output_data='~/ResearchMIT/mixtures/sound-examples';
LF=40;HF=2000;NFB=30;
DUR=500/1000; %100ms
FS0=1000;
IS_SHOW=true;

disp('reading corpus...');
cd (base_corpus);
files=dir('s*.wav');

MYRMS=10; % rms of the mixture!!! important
%%

NF=length(files);
DATA=cell(ITER,1);
cnt=0;
for mm=1:length(Ms),
    M=Ms(mm);
    
    
    for KK=1:ITER,
        
        if mod(KK,30)==1
            fprintf('\tnow %d of %d (%g %%)\n',KK,ITER,KK/ITER*100);
        end
        mixme=cell(M,1);
        myCgrm=cell(M,1);
        mynames=cell(M,1);
        cd (base_corpus);
        mlength=inf;
        for I=1:M,
            iD=randi(NF,1,1);
            
            fname=files(iD).name;
            info=audioinfo(fname);
            smpls=info.TotalSamples;
            fs=info.SampleRate;
            assert(smpls-fs*DUR>0);
            mypos=randi(smpls-fs*DUR,1,1);
            myrange=[mypos, mypos+fs*DUR];
            
            [Y, FS]=audioread(fname,myrange);
            
            %assert(FS==44100);
            Y=double(Y);
            Y=0.03*(Y/sqrt(mean((Y.^2))));
            R=resample(Y,FS0,FS);
            lengthR=length(R)/FS0;
            mixme{I}=0.03*(R/sqrt(mean((R.^2))));
            mlength=min(mlength,lengthR);
            mynames{I}=fname;
        end
        FS=FS0;
        
        soundM=zeros(round(FS0*mlength),1);
        for I=1:M,
            mixme{I}=mixme{I}(1:round(FS0*mlength));
            mixme{I}=mixme{I}*(10^(-MYRMS*(I-1)/20)); %note we make sounds uneven!
            
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
        maskM=(10^(-MYRMS/20))*(myCgrm{1}.^2)>(myCgrm{2}.^2);
        Y1=collapse_subbands(myCgrmM.*(maskM==0),FilterBank);
        Y2=collapse_subbands(myCgrmM.*(maskM==1),FilterBank);
        Y1=0.03*(Y1/sqrt(mean((Y1.^2))));
        Y2=0.03*(Y2/sqrt(mean((Y2.^2))));
        
        
        cnt=cnt+1;
        DATA{cnt}.y1=myCgrm{1};
        DATA{cnt}.y2=myCgrm{2};
        DATA{cnt}.ym=myCgrmM;
        
        DATA{cnt}.mask=maskM;
        DATA{cnt}.mixture=-log2(abs(myCgrmM.^2+eps));
        
        if cnt==1
            DATA{cnt}.DUR=DUR;
            DATA{cnt}.SubbandFrequencies=SubbandFrequencies;
            DATA{cnt}.time=(1:size(myCgrmM,1))/FS;
            DATA{cnt}.timeR=(1:length(Y1))/FS0;
            
            DATA{cnt}.num_featuresX=numel(myCgrmM);
            DATA{cnt}.num_featuresY=numel(myCgrmM);
        end
        assert(numel(DATA{cnt}.y1)==numel(DATA{cnt}.y2));assert(numel(DATA{cnt}.y1)==numel(DATA{cnt}.mask));assert(numel(DATA{cnt}.mask)==numel(DATA{cnt}.mixture));
        
        if IS_SHOW
            figure(1);
            subplot(3,2,1);imagesc(log(myCgrmM.^2)');axis xy;title('mixture');
            subplot(3,2,2);imagesc(maskM');axis xy; title('mask');
            subplot(3,2,3);imagesc(log(myCgrm{1}.^2)');axis xy; title('original 1');
            subplot(3,2,4);imagesc(log(myCgrm{2}.^2)');axis xy; title('original 2');
            
            subplot(3,2,5);imagesc(log((myCgrm{1}.*(maskM==0)).^2 )');axis xy;title('masked mixture 1');
            subplot(3,2,6);imagesc(log((myCgrm{2}.*(maskM==1)).^2 )');axis xy;title('masked mixture 2');
            
            p = audioplayer(soundM, FS);p.play
            pause(1.5*length(soundM)/FS +0.5);
            %          p = audioplayer(mixme{1}, FS);p.play
            %          pause(1.1*length(mixme{1})/FS);
            %          p = audioplayer(mixme{2}, FS);p.play
            %          pause(1.1*length(mixme{2})/FS);
            
            
            
            p1 = audioplayer(Y1, FS);p1.play
            pause(1.1*length(Y1)/FS +0.5);
            p2 = audioplayer(Y2, FS);p2.play
            pause(1.1*length(Y2)/FS+ 0.5);
            %        cd (output_data);
            %
            %         audiowrite('mixture.wav',soundM,FS);
            %         audiowrite('seperated-1.wav',Y1,FS);
            %         audiowrite('seperated-2.wav',Y2,FS);
            pause();
            
        end
        
        
    end
end
%%
disp('reformating input matrices ...');
SubbandFrequencies=DATA{1}.SubbandFrequencies;
time=DATA{1}.time;
timeR=DATA{1}.time;
num_featuresX=DATA{1}.num_featuresX;
num_featuresY=DATA{1}.num_featuresY;
X=nan(ITER,num_featuresX);
Y=nan(ITER,num_featuresY);
for  I=1:ITER,
    if mod(I,30)==1
        fprintf('\tnow %d of %d (%g %%)\n',I,ITER,I/ITER*100);
    end
    X(I,:)=DATA{I}.mixture(:);
    posp=DATA{I}.mask(:)==1;
    posn=DATA{I}.mask(:)~=1;
    Y(I,posp)=1;
    Y(I,posn)=-1;
end

[Xn,MuX,SigX]=zscore(X);
%Xn=X;
Yn=Y; %allready normalized
save('temp.mat');
%%
% ********** Examples of using CCA *********
options.PrjX = 1;
options.PrjY = 1;
options.RegX = .7;
options.RegY = .7;
[W_x, W_y, corr_list] = CCA(Xn', Yn', options);
X_p = Xn* W_x;
Y_p = Yn* W_y;
%% 
%%
Ap=16;
figure(100);clf;
   
for J=1:min(size(W_x,2),Ap),
    subplot(4,4,J);
    %plot(sum(reshape(((W_x(:,J))),[length(timeR) length(SubbandFrequencies)])',2))
     imagesc(reshape(((W_x(:,J))),[length(timeR) length(SubbandFrequencies)])');
     
   % title(r(J));
end
figure(101);clf;
for J=1:min(size(W_y,2),Ap),
    subplot(4,4,J);
    imagesc(reshape(((W_y(:,J))),[length(timeR) length(SubbandFrequencies)])');
    %plot(sum(reshape(((W_y(:,J))),[length(timeR) length(SubbandFrequencies)])',2))
    %plot(sum(reshape(((W_y(:,J))),[length(timeR) length(SubbandFrequencies)])',2))
    axis off;
   % title(r(J));
end

figure(100);
%%

[A B r U V] = canoncorr(Xn, Yn);
%%
Ap=100;
figure(100);clf;
   
for J=1:min(size(A,2),Ap),
    subplot(10,10,J);
    imagesc(reshape(A(:,J)',[length(timeR) length(SubbandFrequencies)])')
    axis off;
   % title(r(J));
end
figure(101);clf;
for J=1:min(size(B,2),Ap),
    subplot(10,10,J);
    imagesc(reshape(B(:,J)',[length(timeR) length(SubbandFrequencies)])')
    axis off;
   % title(r(J));
end
%%
addpath


%%
figure(3);clf;
% ktypeX='gaussian';paramX=100;
  ktypeX='polynomial';paramX=1;

display('calculating kernels (X)...');
KX=TOOLk_calc_kernel(Xn',Xn',ktypeX,paramX)';
subplot(2,2,1);imagesc(KX-diag(diag(KX)));

  ktypeY='polynomial';paramY=1;
%ktypeY='gaussian';paramY=100;
display('calculating kernels (Y)...');
KY=TOOLk_calc_kernel(Y',Y',ktypeY,paramY)';
subplot(2,2,2);imagesc(KY-diag(diag(KY)));

%%
% 
% ktypeY='polynomial';paramY=1;
% display('calculating kernels (Y)...');
% KY=TOOLk_calc_kernel(Y',Y',ktype,paramY)';
% 


kappa=0.1;
eta=1;

display('doing kCCA...');
[nalpha, nbeta, rk, KXC, KYC] = TOOLk_kcanonca_reg_ver2(KX,KY,kappa ,eta);

IBr=  mean(   (    KXC*(nalpha(:,end:-1:1)) )  .^2  )' ;
IBlambda=1-rk(end:-1:1).*rk(end:-1:1);
IBbetaC=1./(1-IBlambda);
IBcurveBasis=TOOLk_calc_IB_original(rk(end:-1:1),IBbetaC) ;
figure(4);clf;plot(IBcurveBasis(:,1),IBcurveBasis(:,2),'o-');

%%
% IBbeta=IBbetaC;
%
% [weights,IBX,IBY]=TOOLk_calc_weights_IBX_IBY(KXC,nalpha(:,end:-1:1),rk(end:-1:1),IBbeta);
%     IBcurve=[IBcurve;[IBX,IBY]];
%     if (sum(weights>0)<=1) && (abs(weights(1))<1e-4)
%         sprintf('Set weights from %g to 0', weights(1));
%         weights(1)=0;
%
%     end
%






