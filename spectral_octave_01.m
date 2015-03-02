% compute features
%cd 'C:\Users\user\Dropbox (PPCA)\Research MIT\mixtures'
close all;clc;clear all;
NFB=30;

Ms=[2];
MN=length(Ms);
ITER=1000;

LF=40;HF=8000;

figure(1);clf;



myhist_all=zeros(NFB,NFB);

%mfname='data-inst/mix-v1.M.2.3.wav';
mfname='~/ResearchMIT/assorted/Octave/OI3-02 LONG SEQUENCE-nochange.wav';

display(sprintf('%s',mfname));

[tsS,fs]=audioread(mfname, [1+0*44000 ,4*44000]);
%ts=tsS(:,1)+tsS(:,2);
%p = audioplayer(0.1*tsS(:,1), fs);p.play
%ts=tsS(:,1);

%ts=ts/sqrt(mean(ts.*ts));
[FilterBank,SubbandFrequencies]=make_erb_cos_filters(size(tsS,1),fs,NFB-2,LF,HF);

CgrmS=cell(3,1);
CgrmS{1}= generate_subbands(tsS(:,1)/mean(sqrt(tsS(:,1).^2)),FilterBank);
CgrmS{2}= generate_subbands(tsS(:,2)/mean(sqrt(tsS(:,2).^2)),FilterBank);
CgrmS{3}= generate_subbands((tsS(:,1)+tsS(:,2))/mean(sqrt((tsS(:,1)+tsS(:,2)).^2)),FilterBank);

CgrmA=[CgrmS{1},CgrmS{2},CgrmS{3}];

myCgrmS=cell(3,1);
for I=1:3,
    %myCgrmS{I}=sqrt((CgrmS{I}.*CgrmS{I}));
    myCgrmS{I}=abs(hilbert(CgrmS{I}));
end

myCgrm=[myCgrmS{1},myCgrmS{2},myCgrmS{3}];

NFB3=3*NFB;
myhist=zeros(NFB3,NFB3);
JUMP= 0;
for I=1:NFB3,
    fprintf('I=%d of %d\n',I,NFB3);
    for J=1:NFB3,
        
        myhist(I,J)=sum( (myCgrm(:,I)).*( myCgrm(:,J)))/length(myCgrm(:,J));
        
    end
end

myhist=myhist/sum(sum(myhist));
% marg=sum(myhist);
% myhist_marg=repmat(marg,size(myhist,1),1).*repmat(marg',1,size(myhist,1));
% pmi=log2(myhist+eps)./(log2(myhist_marg+eps));
%%
subplot(2,3,1);
%figure(2);
colormap (bone);
imagesc(log2((myhist+0.00000)+eps));
title('coherence');axis xy;hold on;
colorbar;
plot([1 NFB3],[NFB,NFB],'y-','LineWidth',3);
plot([NFB,NFB],[1 NFB3],'y-','LineWidth',3);
plot([1 NFB3],[2*NFB,2*NFB],'y-','LineWidth',3);
plot([2*NFB,2*NFB],[1 NFB3],'y-','LineWidth',3);
plot([1 NFB3],[NFB,NFB],'y-','LineWidth',3);
plot([NFB,NFB],[1 NFB3],'y-','LineWidth',3);
plot([1 NFB3],[2*NFB,2*NFB],'y-','LineWidth',3);
plot([2*NFB,2*NFB],[1 NFB3],'y-','LineWidth',3);
% plot([1 NFB3],14+[0,0],'b--','LineWidth',1);
% plot(14+[0,0],[1 NFB3],'b--','LineWidth',1);
% plot([1 NFB3],22+[0,0],'r--','LineWidth',1);
% plot(22+[0,0],[1 NFB3],'r--','LineWidth',1);
set(gca,'Xtick',[14,22,14+NFB,22+NFB,14+NFB*2,22+NFB*2]);
set(gca,'XtickLabel',{'LL','LH','RL','RH','BL','BH'});
set(gca,'Ytick',[14,22,14+NFB,22+NFB,14+NFB*2,22+NFB*2]);
set(gca,'YtickLabel',{'LL','LH','RL','RH','BL','BH'});


%%
subplot(2,3,3);


colormap(bone);



imagesc(nm);title('nm');axis xy;hold on;

plot([1 NFB3],[NFB,NFB],'y-','LineWidth',3);
plot([NFB,NFB],[1 NFB3],'y-','LineWidth',3);
plot([1 NFB3],[2*NFB,2*NFB],'y-','LineWidth',3);
plot([2*NFB,2*NFB],[1 NFB3],'y-','LineWidth',3);
plot([1 NFB3],[NFB,NFB],'y-','LineWidth',3);
plot([NFB,NFB],[1 NFB3],'y-','LineWidth',3);
plot([1 NFB3],[2*NFB,2*NFB],'y-','LineWidth',3);
plot([2*NFB,2*NFB],[1 NFB3],'y-','LineWidth',3);
set(gca,'Xtick',[14,22,14+NFB,22+NFB,14+NFB*2,22+NFB*2]);
NAMS={'LL','LH','RL','RH','BL','BH'};
NAMS={'Left-Low','Left-High','Right-Low','Right-High','Binaural-Low','Binaural-high'};
set(gca,'XtickLabel',NAMS);
set(gca,'Ytick',[14,22,14+NFB,22+NFB,14+NFB*2,22+NFB*2]);
set(gca,'YtickLabel',NAMS);


colorbar;

%%

%%
figure(1);


[C, ~, ~] = SpectralClustering(nm, 2, 2);



C=full(C);
%tpos=(sum(nm)>max(sum(nm))/PTH)';
pos1=(C(:,1)==1); %& tpos ;
pos2=(C(:,2)==1); %& tpos ;
pos3=[find(pos1);find(pos2)];

Cgrm1=zeros(size(myCgrm));Cgrm1(:,pos1)=myCgrm(:,pos1);
Cgrm2=zeros(size(myCgrm));Cgrm2(:,pos2)=myCgrm(:,pos2);

subplot(2,3,2);
%imagesc(log2(myhist(pos3,pos3)+eps));title('re-ordered hist');axis xy;hold on;
imagesc((nm(pos3,pos3)));title('re-ordered hist');axis xy;hold on;
plot([sum(pos1),sum(pos1)],[0,length(pos3)],'y-','LineWidth',3);
plot([0,length(pos3)],[sum(pos1),sum(pos1)],'y-','LineWidth',3);


colorbar;

%%
subplot(2,3,4);
imagesc(log2(abs((myCgrm'+eps))));title('spectrogram');axis xy;hold on;
plot([0,size(myCgrm,1)],[NFB,NFB],'y-','LineWidth',3);
plot([0,size(myCgrm,1)],[2*NFB,2*NFB],'y-','LineWidth',3);
whr=[14, 22,14+NFB,22+NFB,14+NFB*2,22+NFB*2];
whc=pos1(whr);
CLRS={'b-','r-'};

for I=1:length(whr)
    plot([0,size(myCgrm,1)],[whr(I)-2,whr(I)-2],CLRS{whc(I)+1},'LineWidth',2);
    plot([0,size(myCgrm,1)],[whr(I)+2,whr(I)+2],CLRS{whc(I)+1},'LineWidth',2);
    
end
set(gca,'Ytick',[14,22,14+NFB,22+NFB,14+NFB*2,22+NFB*2]);
set(gca,'YtickLabel',NAMS);

% plot([0,size(myCgrm,1)],[14 14],'b-','LineWidth',2);
% plot([0,size(myCgrm,1)],[22 22],'r-','LineWidth',2);
% plot([0,size(myCgrm,1)],14+[NFB,NFB],'b-','LineWidth',2);
% plot([0,size(myCgrm,1)],14+[2*NFB,2*NFB],'b-','LineWidth',2);
% plot([0,size(myCgrm,1)],22+[NFB,NFB],'r-','LineWidth',2);
% plot([0,size(myCgrm,1)],22+[2*NFB,2*NFB],'r-','LineWidth',2);
%
%
%%
subplot(2,3,3);
whr=[14, 22,14+NFB,22+NFB,14+NFB*2,22+NFB*2];
whc=pos1(whr);
CLRS={'b-','r-'};

for I=1:length(whr)
    plot([0,NFB3],[whr(I)-2,whr(I)-2],CLRS{whc(I)+1},'LineWidth',2);
    plot([0,NFB3],[whr(I)+2,whr(I)+2],CLRS{whc(I)+1},'LineWidth',2);
    plot([whr(I)-2,whr(I)-2],[0,NFB3],CLRS{whc(I)+1},'LineWidth',2);
    plot([whr(I)+2,whr(I)+2],[0,NFB3],CLRS{whc(I)+1},'LineWidth',2);
end
set(gca,'Ytick',[14,22,14+NFB,22+NFB,14+NFB*2,22+NFB*2]);
set(gca,'YtickLabel',NAMS);





%%

Cgrm1p=Cgrm1+0*0.00051*Cgrm2;
Cgrm2p=Cgrm2+0*0.00051*Cgrm1;
subplot(2,3,5);
nm1p=zeros(size(nm));nm1p(pos1,pos1)=((nm(pos1,pos1)));
nm2p=zeros(size(nm));nm2p(pos2,pos2)=((nm(pos2,pos2)));
hold off;plot(sum(nm1p),'y-','LineWidth',3);title('spectrogram1');hold on;
tmax=max(sum(nm1p));
plot([14 14],[0 tmax],'b:');
plot([22 22],[0 tmax],'r:');
plot(NFB+[14 14],[0 tmax],'b:');
plot(NFB+[22 22],[0 tmax],'r:');
plot(2*NFB+[14 14],[0 tmax],'b:');
plot(2*NFB+[22 22],[0 tmax],'r:');
plot(NFB+[0 0],[0 tmax],'g--','LineWidth',2);
plot(2*NFB+[0 0],[0 tmax],'g--','LineWidth',2);

subplot(2,3,6);
hold off;plot(sum(nm2p));title('spectrogram2');hold on;
plot(sum(nm2p),'y-','LineWidth',3);title('spectrogram2');hold on;
tmax=max(sum(nm2p));
plot([14 14],[0 tmax],'b:');
plot([22 22],[0 tmax],'r:');
plot(NFB+[14 14],[0 tmax],'b:');
plot(NFB+[22 22],[0 tmax],'r:');
plot(2*NFB+[14 14],[0 tmax],'b:');
plot(2*NFB+[22 22],[0 tmax],'r:');
plot(NFB+[0 0],[0 tmax],'g--','LineWidth',2);
plot(2*NFB+[0 0],[0 tmax],'g--','LineWidth',2);


%%

% subplot(2,3,5);
% imagesc(log2(sqrt((Cgrm1p.*Cgrm1p))'+eps));title('spectrogram1');axis xy;hold on;
% plot([0,size(myCgrm,1)],[NFB,NFB],'y-','LineWidth',3);
% plot([0,size(myCgrm,1)],[2*NFB,2*NFB],'y-','LineWidth',3);
% plot([0,size(myCgrm,1)],[14 14],'b-','LineWidth',2);
% plot([0,size(myCgrm,1)],[22 22],'r-','LineWidth',2);
% plot([0,size(myCgrm,1)],14+[NFB,NFB],'b-','LineWidth',2);
% plot([0,size(myCgrm,1)],14+[2*NFB,2*NFB],'b-','LineWidth',2);
% plot([0,size(myCgrm,1)],22+[NFB,NFB],'r-','LineWidth',2);
% plot([0,size(myCgrm,1)],22+[2*NFB,2*NFB],'r-','LineWidth',2);
% subplot(2,3,6);
% imagesc(log2(sqrt((Cgrm2p.*Cgrm2p))'+eps));title('spectrogram2');axis xy;hold on;
%
% plot([0,size(myCgrm,1)],[NFB,NFB],'y-','LineWidth',3);
% plot([0,size(myCgrm,1)],[2*NFB,2*NFB],'y-','LineWidth',3);
% plot([0,size(myCgrm,1)],[14 14],'b-','LineWidth',2);
% plot([0,size(myCgrm,1)],[22 22],'r-','LineWidth',2);
% plot([0,size(myCgrm,1)],14+[NFB,NFB],'b-','LineWidth',2);
% plot([0,size(myCgrm,1)],14+[2*NFB,2*NFB],'b-','LineWidth',2);
% plot([0,size(myCgrm,1)],22+[NFB,NFB],'r-','LineWidth',2);
% plot([0,size(myCgrm,1)],22+[2*NFB,2*NFB],'r-','LineWidth',2);

drawnow;
Y1S=cell(3,1);
Y2S=cell(3,1);

CgrmA1=zeros(size(CgrmA));CgrmA1(:,pos1)=CgrmA(:,pos1);
CgrmA2=zeros(size(CgrmA));CgrmA2(:,pos2)=CgrmA(:,pos2);




for I=1:3,
    
    Y1S{I}=0.5*collapse_subbands(CgrmA1(:,(1+(I-1)*NFB):(I*NFB)),FilterBank);
    Y2S{I}=0.5*collapse_subbands(CgrmA2(:,(1+(I-1)*NFB):(I*NFB)),FilterBank);
end

p = audioplayer(Y1S{3}/max(Y1S{3}), fs);p.play
pause(1.5*length(Y1S{3})/fs);

p = audioplayer(Y2S{3}/max(Y2S{3}), fs);p.play
pause(1.5*length(Y2S{3})/fs);

%p = audioplayer([Y2S{1}/max(Y2S{1}),Y2S{2}/max(Y2S{2})], fs);p.play
%pause(1.5*length(Y2S{I})/fs);

% for I=1:3,
%     p = audioplayer(Y1S{I}/max(Y1S{I}), fs);p.play
%     pause(1.5*length(Y1S{I})/fs);
%
%     p = audioplayer(Y2S{I}/max(Y2S{I}), fs);p.play
%     pause(1.5*length(Y2S{I})/fs);
%
% end

% for tt=1:2,
%
%     p = audioplayer(0.1*ts, fs);p.play
%     pause(1.5*length(ts)/fs);
%     p1 = audioplayer(0.1*Y1, fs);p1.play
%     pause(1.1*length(ts)/fs);
%     p2 = audioplayer(0.1*Y2, fs);p2.play
%     pause(1.1*length(ts)/fs);
%     pause
%
%     audiowrite('SPEC-mix.wav',0.1*ts,fs);
%     audiowrite('SPEC-unmix1.wav',0.1*Y1,fs);
%     audiowrite('SPEC-unmix2.wav',0.1*Y2,fs);
% end
%
% figure(2);
% subplot(3,3,m);
% imagesc(SubbandFrequencies,SubbandFrequencies,myhist_all);
% title(sprintf('num spks %d',M));
