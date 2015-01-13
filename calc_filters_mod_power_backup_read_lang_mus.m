% assert(1==0);
clear all;close all;clc;
disp(sprintf('reading and formatting data...'));
% load('FEATURES-N-v4.mat');
% load('FEATURES-NP-v2.mat'); %timit mixture speech
%  load('FEATURES-cello-v1.mat');
% load('FEATURES-AUDIO-BOOKS-ALL');
% load('FEATURES-AUDIO-BOOKS-ENGLISH.mat');
% load('FEATURES-AUDIO-BOOKS-ITAL.mat');
%  load('FEATURES-AUDIO-BOOKS-TOG.mat');
% load('FEATURES-AUDIO-BOOKS-HUN.mat');
load('FEATURES-CELLO2.mat');

Nab=2;
FEATURESab=cell(Nab,1);
A=1;B=2;
FEATURESab{1}=FEATURES{A};
FEATURESab{2}=FEATURES{B};
FEATURESab{1}(isnan(FEATURESab{1}))=0;
FEATURESab{2}(isnan(FEATURESab{2}))=0;

INFOab=cell(Nab,1);
INFOab{1}=INFO{A};
INFOab{2}=INFO{B};


% addpath('C:\Users\user\Dropbox (PPCA)\Research MIT\mixtures\libsvm');
addpath('C:\Users\user\Dropbox (PPCA)\Research MIT\mixtures\libplsda');

NTRAIN=round(size(FEATURESab{1},1)*9/10);

vals=[];labels=[];
for I=1:Nab
    vals=[vals;FEATURESab{I}(1:NTRAIN,:)];
    labels=[labels;ones(size(FEATURESab{I}(1:NTRAIN,:),1),1)*I];
    
end

info=cell(length(labels),1);cnt=0;
for I=1:Nab
    for J=1:NTRAIN
        cnt=cnt+1;
        info{cnt,1}.fname=INFOab{I}{J}.fname;
        info{cnt,1}.range=INFOab{I}{J}.range;
        %         info{cnt,1}.Hz_mod_cfreqs=INFOab{I}{J,1}.Hz_mod_cfreqs;
        %         info{cnt,1}.audio_cutoffs_Hz=INFOab{I}{J,1}.audio_cutoffs_Hz;
        Hz_mod_cfreqs=INFOab{I}{J,1}.Hz_mod_cfreqs;
        audio_cutoffs_Hz=INFOab{I}{J,1}.audio_cutoffs_Hz;
        fs=INFOab{I}{J,1}.fs;
        info{cnt,1}.audio=INFOab{I}{J}.audio;
        info{cnt,1}.fs=INFOab{I}{J}.fs;
        
    end
end
assert(length(labels)==length(info));
%%
[vals,MU,SIGMA] = zscore(vals);


valsT=[];labelsT=[];
for I=1:Nab
    valsT=[valsT;FEATURESab{I}((NTRAIN+1):end,:)];
    labelsT=[labelsT;ones(size(FEATURESab{I}((NTRAIN+1):end,:),1),1)*I];
end
valsT=(valsT-repmat(MU,size(valsT,1),1))./repmat(SIGMA,size(valsT,1),1);
valsT(isnan(valsT))=0;

infoT=cell(length(labelsT),1);cnt=0;
for I=1:Nab
    for J=(NTRAIN+1):(size(FEATURESab{I},1))
        cnt=cnt+1;
        infoT{cnt,1}.fname=INFOab{I}{J}.fname;
        infoT{cnt,1}.range=INFOab{I}{J}.range;
        %         infoT{cnt,1}.Hz_mod_cfreqs=INFOab{I}{J,1}.Hz_mod_cfreqs;
        %         infoT{cnt,1}.audio_cutoffs_Hz=INFOab{I}{J,1}.audio_cutoffs_Hz;
        Hz_mod_cfreqs=INFOab{I}{J,1}.Hz_mod_cfreqs;
        audio_cutoffs_Hz=INFOab{I}{J,1}.audio_cutoffs_Hz;
        
    end
end
assert(length(labelsT)==length(infoT));
%%
Xp=vals;
yp=labels*2-3;
Ap=20; %maximal number of dimensions in projections
XpT=valsT;
ypT=labelsT*2-3;

disp(sprintf('doing cross-validation...'));
CV=plscv(Xp,yp,Ap);

disp(sprintf('doing pls pn all data...'));
PLS=pls(Xp,yp,Ap);
acc=sum(sign(CV.Ypred(:,:))==repmat(yp,1,size(CV.Ypred,2)))/length(yp)
disp(sprintf('predict on test set...'));
[ypredT]=plsval(PLS,XpT,ypT,Ap);
accT=sum(sign(ypredT)==ypT)/length(ypT);

%%
disp(sprintf('displaying data projected to sub spaces...'));
figure(55);clf;
NTS=5;
for I=1:NTS,
    for J=(I):NTS,
        subplot(NTS,NTS,(I-1)*NTS+J);
        if I==J
            plot(PLS.X_scores(yp>0,J)+randn(size(PLS.X_scores(yp>0,J))),PLS.X_scores(yp>0,I)+randn(size(PLS.X_scores(yp>0,I))),'.b');hold on;
            plot(PLS.X_scores(yp<0,J)+randn(size(PLS.X_scores(yp<0,J))),PLS.X_scores(yp<0,I)+randn(size(PLS.X_scores(yp<0,I))),'.r');
        else
            plot(PLS.X_scores(yp>0,J),PLS.X_scores(yp>0,I),'.b');hold on;
            plot(PLS.X_scores(yp<0,J),PLS.X_scores(yp<0,I),'.r');
        end
        xlabel(sprintf('%d',J));
        ylabel(sprintf('%d',I));
    end
end

%%

disp(sprintf('showing filters...'));

cmp=zeros(64,3)
for I=1:32,
    cmp(32+I,:)=[0,0,I/32];
    cmp(I,:)=[1-I/32,0,0];
    
end
figure(17);clf;figure(18);clf;
for I=1:min(15,Ap),
    figure(17);subplot(4,4,I);
    wf=PLS.X_loadings(:,I);
    
    
    xlgnd=audio_cutoffs_Hz;
    ylgnd=Hz_mod_cfreqs;
    wfr=reshape(wf,[length(audio_cutoffs_Hz),length(Hz_mod_cfreqs)]);
    
    imagesc(ylgnd,xlgnd,wfr);axis xy;
    colormap(cmp)
    
    ylabel('Fq (Hz)');
    xlabel('mod Fq(Hz)');
    title(sprintf('%d',I));
    
end

figure(17);subplot(4,4,16);

plot(acc*100,'kx-','LineWidth',2);hold on;
plot([1 Ap],[accT,accT]*100,'g--','LineWidth',2);hold on;
% plot([15 15],[0,100],'r--','LineWidth',1);hold on;
axis([1 Ap 50 100]);
title('Labeling accuracy');
xlabel('Number of coefficents');

figure(18);
ylabel('Accuracy (%%)');
h=legend('cv-train','all-test','Location','NorthEastOutside');
set(h,'FontSize',8);


plot(acc*100,'kx-','LineWidth',2);hold on;
plot([1 Ap],[accT,accT]*100,'g--','LineWidth',2);hold on;
% plot([15 15],[0,100],'r--','LineWidth',1);hold on;
axis([1 Ap 50 100]);
title('Labeling accuracy');
xlabel('Number of coefficents');


ylabel('Accuracy (%%)');
h=legend('cv-train','all-test','Location','SouthEast');
set(h,'FontSize',8);
%%
disp(sprintf('displaying raw data as a matrix...'));

figure(1);clf;imagesc([FEATURES{2};FEATURES{1}]);axis xy
%%

WHO=[1,2];
for J=1:length(WHO),
    figure(30+J-1);clf;
    MAXSOUNDS=25;
    disp(sprintf('displaying examples of stimuli'));
    for I=1:MAXSOUNDS,
        subplot(floor(sqrt(MAXSOUNDS)+0.5),floor(sqrt(MAXSOUNDS)+0.5),I)
        idx=find(labels==WHO(J));
        xlgnd=audio_cutoffs_Hz;
        ylgnd=Hz_mod_cfreqs;
        
        imagesc(ylgnd,xlgnd, (reshape(vals(idx(I),:),[length(xlgnd),length(ylgnd)]))');axis xy;
    end
end

%%

figure(2);clf;

hold on;I1=2;I2=1;
pos=(labels~=1);plot(PLS.X_scores(pos,I1),PLS.X_scores(pos,I2),'.r');
pos=(labels==1);plot(PLS.X_scores(pos,I1),PLS.X_scores(pos,I2),'.b');

[val,idx]=sort(PLS.X_scores(:,I1));



%%
disp(sprintf('concatenating typical audio examples...'));

figure(2);clf;

hold on;
WHO=[1:Ap]
% WHO=[4]
MAXSOUNDS=25;
for I1=WHO,
    % I1=6;
    % I2=1;
    I2=min(Ap,I1+1);
    % pos=(labels~=1);plot(PLS.X_scores(pos,I1),PLS.X_scores(pos,I2),'.r');
    % pos=(labels==1);plot(PLS.X_scores(pos,I1),PLS.X_scores(pos,I2),'.b');
    
    [val,idx]=sort(PLS.X_scores(:,I1));
    Nplay=400;
    RANGES={[1:Nplay],[length(idx):-1:(length(idx)-Nplay)]};
    
    for k=1:length(RANGES),
        %     if k==1
        %         playme=cos(2*pi*(1:8000)/fs*300)';
        %     else
        %         playme=0.3*cos(2*pi*(1:8000)/fs*300)'+0.3*cos(2*pi*(1:8000)/fs*450)';
        %     end
        playme=zeros(fs,1);
        
        figure(100+k);clf
        cnt=1;
        hashy=inf;
        MYHASHES=[];
        for J=RANGES{k},
            if cnt>MAXSOUNDS;
                break
            end
            if sum(vals(idx(J),:))~=hashy
                subplot(floor(sqrt(MAXSOUNDS)+0.5),floor(sqrt(MAXSOUNDS)+0.5),min(cnt,MAXSOUNDS)); cnt=cnt+1;
                hashy= sum(vals(idx(J),:));
                xlgnd=audio_cutoffs_Hz;
                
                ylgnd=Hz_mod_cfreqs;
                imagesc(xlgnd,ylgnd,reshape(vals(idx(J),:),[length(xlgnd),length(ylgnd)]));axis xy;
                axis xy;
                axis off;
                myrng=info{idx(J)}.range;
                %                 myrng=[round(max(info{idx(J)}.range)*1/3),round(max(info{idx(J)}.range)*2/3)];

                % do this if these are not saved
%                 [Y, fs]=audioread(info{idx(J)}.fname,myrng );

                Y=info{idx(J)}.audio;
                fs=info{idx(J)}.fs;
                
                %                 HASH=sum(Y)+sum(round(Y*1000))/1000;
                %                 if sum(find(MYHASHES==HASH))==0
                %                     MYHASHES=[MYHASHES,HASH];
                
                Yh=Y.*hanning(length(Y));
                n3=round(length(Yh)/3);
                playmeold=playme;
                playme=zeros(length(playmeold)+3*n3,1);
                playme(1:length(playmeold))=playmeold;
                playme((length(playmeold)-n3+1):(length(playmeold)-n3+length(Yh)))=playme((length(playmeold)-n3+1):(length(playmeold)-n3+length(Yh)))+Yh;
                %                 else
                %                     fprintf('ignoring %g\n',HASH);
                %
                %                 end
            end
        end
        
        
        drawnow;
        p=audioplayer(playme,fs);p.play;
        ofname=sprintf('NPLS-concat-d-%d.direc-%d.wav',I1,k);
        
        audiowrite(ofname,playme,fs)
        pause(0.8*length(playme)/fs);
        pause
    end
end
%%
