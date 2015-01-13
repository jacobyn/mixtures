% assert(1==0);
clear all;close all;clc;
disp(sprintf('reading and formatting data...'));
% load('FEATURES-N-v4.mat'); 
load('FEATURES-NI-v1.mat');
Nab=2;
FEATURESab=cell(Nab,1);
A=1;B=2;
FEATURESab{1}=FEATURES{A};
FEATURESab{2}=FEATURES{B};
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

infoT=cell(length(labelsT),1);cnt=0;
for I=1:Nab
    for J=(NTRAIN+1):(size(FEATURESab{I},1))
        cnt=cnt+1;
        infoT{cnt,1}.fname=INFOab{I}{J}.fname;
        infoT{cnt,1}.range=INFOab{I}{J}.range;
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
    xlgnd=1000*(1:JMP:(WIND*FS0))./FS0;
    ylgnd=SubbandFrequencies;
    wfr=reshape(wf,[NFB,NTP]);
    %     figure(11);clf;
    %      subplot(2,2,1);
    imagesc(xlgnd,ylgnd,wfr);hold on;axis xy;
    
    
    colormap(cmp)
    
    ylabel('Fq (Hz)');
    xlabel('time(ms)');
    title(sprintf('%d',I));
    
    figure(18);subplot(4,4,I);
    plot(ylgnd,sum(reshape(PLS.X_loadings(:,I),[NFB,NTP]),2));
    xlabel('Fq (Hz)');
    %     ylabel('weight');
    title(sprintf('marginal of comp %d',I));
    
end

figure(17);subplot(4,4,16);

plot(acc*100,'kx-','LineWidth',2);hold on;
plot([1 Ap],[accT,accT]*100,'g--','LineWidth',2);hold on;
% plot([15 15],[0,100],'r--','LineWidth',1);hold on;
axis([1 Ap 50 100]);
title('Labeling accuracy');
xlabel('Number of coefficents');


ylabel('Accuracy (%%)');
h=legend('cv-train','all-test','Location','NorthEastOutside');
set(h,'FontSize',8);
%%
disp(sprintf('displaying raw data as a matrix...'));

figure(1);clf;imagesc([FEATURES{2};FEATURES{1}]);axis xy

figure(17);

%%
% [x,fs]=wavread(INFO{1}{1}.fname); p = audioplayer(x([min(INFO{1}{1}.range):(max(INFO{1}{1}.range)),min(INFO{1}{1}.range):(max(INFO{1}{1}.range))]), fs);p.play



%%

figure(2);clf;

hold on;I1=2;I2=1;
pos=(labels~=1);plot(PLS.X_scores(pos,I1),PLS.X_scores(pos,I2),'.r');
pos=(labels==1);plot(PLS.X_scores(pos,I1),PLS.X_scores(pos,I2),'.b');

[val,idx]=sort(PLS.X_scores(:,I1));


%
% %%
% playme=[];
% for J=1:Nplay,
%     [Y, fs]=audioread(info{idx(J)}.fname, info{idx(J)}.range);
%     playme=[playme;Y];
% end
% p=audioplayer(playme,fs);p.play;
%
% %%
% playme=[];
% for J=(length(idx)-Nplay):length(idx),
%     [Y, fs]=audioread(info{idx(J)}.fname, info{idx(J)}.range);
%     playme=[playme;Y];
% end
% p=audioplayer(playme,fs);p.play;

%%
disp(sprintf('concatenating typical audio examples...'));

figure(2);clf;

hold on;
WHO=[1:Ap]
% WHO=[4]
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
    playme=zeros(4000,1);
    
    figure(100+k);clf
    cnt=1;
    hashy=inf;
    
    for J=RANGES{k},
        if cnt>100;
            break
        end
        if sum(vals(idx(J),:))~=hashy
            subplot(10,10,min(cnt,100)); cnt=cnt+1;
            hashy= sum(vals(idx(J),:));
            xlgnd=1000*(1:JMP:(WIND*FS0))./FS0;
            ylgnd=SubbandFrequencies;
            imagesc(xlgnd,ylgnd,reshape(vals(idx(J),:),[NFB,NTP]));
            axis xy;
            axis off;
            
            
            [Y, fs]=audioread(info{idx(J)}.fname, info{idx(J)}.range);
            Yh=Y.*hanning(length(Y));
            n3=round(length(Yh)/3);
            playmeold=playme;
            playme=zeros(length(playmeold)+3*n3,1);
            playme(1:length(playmeold))=playmeold;
            playme((length(playmeold)-n3+1):(length(playmeold)-n3+length(Yh)))=playme((length(playmeold)-n3+1):(length(playmeold)-n3+length(Yh)))+Yh;
        end
    end
    
    
    drawnow;
    p=audioplayer(playme,fs);p.play;
    ofname=sprintf('concat-d-%d.direc-%d.wav',I1,k);
    audiowrite(ofname,playme,fs)
    pause(0.8*length(playme)/fs);
end
end
%%
disp(sprintf('displaying random points(matrix)...'));

for J=1:2,
    figure(204+J);clf;
    pos=find(labels==J);
    for I=1:100,
        subplot(10,10,I);
        idx=    pos(I);
        [Y, fs]=audioread(info{idx}.fname, info{idx}.range);
        imagesc(xlgnd,ylgnd,reshape(vals(idx,:),[NFB,NTP]));
        axis xy;
        axis off;
        
        
    end
end
%%
disp(sprintf('displaying random points(single)...'));
for J=1:2,
    figure(206+J);clf;
    pos=find(labels==J);
    for I=1:1,
        %     subplot(10,10,I);
        idx=    pos(I);
        [Y, fs]=audioread(info{idx}.fname, info{idx}.range);
        imagesc(xlgnd,ylgnd,reshape(vals(idx,:),[NFB,NTP]));
        axis xy;
        %             axis off;
        ylabel('Fq (Hz)');
        xlabel('time(ms)');
        title(sprintf('num spk= %d',J));
        
    end
end



%%




figure(200);
[Y, fs]=audioread(info{idx(1)}.fname, info{idx(1)}.range);
plot(1000*(1:length(Y))/fs,Y);
%%
figure(201);
imagesc(xlgnd,ylgnd,reshape(vals(idx,:),[NFB,NTP]));axis xy;
ylabel('Fq (Hz)');
xlabel('time(ms)');
%%
figure(210);clf;
cnt=0;
Apm=5;
for I1=1:Apm,
    for I2=1:Apm,
        cnt=cnt+1;
        subplot(Apm,Apm,cnt);
        
        pos=(labels~=1);pos=find(pos);pos=pos(1:10:end);
        plot(PLS.X_scores(pos,I1),PLS.X_scores(pos,I2),'.r');hold on;
        pos=(labels==1);pos=find(pos);pos=pos(1:10:end);
        plot(PLS.X_scores(pos,I1),PLS.X_scores(pos,I2),'.b');hold on;
        
    end
end

