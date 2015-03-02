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

[tsS,fs]=audioread(mfname, [1+0*44000 ,5*44000]);
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
% subplot(2,2,1);
% %figure(2);
% colormap (bone);
% imagesc(log2((myhist+0.00000)+eps));
% title('coherence');axis xy;hold on;
% colorbar;
% plot([1 NFB3],[NFB,NFB],'y-','LineWidth',3);
% plot([NFB,NFB],[1 NFB3],'y-','LineWidth',3);
% plot([1 NFB3],[2*NFB,2*NFB],'y-','LineWidth',3);
% plot([2*NFB,2*NFB],[1 NFB3],'y-','LineWidth',3);
% plot([1 NFB3],[NFB,NFB],'y-','LineWidth',3);
% plot([NFB,NFB],[1 NFB3],'y-','LineWidth',3);
% plot([1 NFB3],[2*NFB,2*NFB],'y-','LineWidth',3);
% plot([2*NFB,2*NFB],[1 NFB3],'y-','LineWidth',3);
%
% set(gca,'Xtick',[14,22,14+NFB,22+NFB,14+NFB*2,22+NFB*2]);
% set(gca,'XtickLabel',{'LL','LH','RL','RH','BL','BH'});
% set(gca,'Ytick',[14,22,14+NFB,22+NFB,14+NFB*2,22+NFB*2]);
% set(gca,'YtickLabel',{'LL','LH','RL','RH','BL','BH'});


%%

%%

%%


nm=(myhist+eps);
[C, ~, ~] = SpectralClustering(nm, 2, 3);



C=full(C);
%tpos=(sum(nm)>max(sum(nm))/PTH)';
pos1=(C(:,1)==1); %& tpos ;
pos2=(C(:,2)==1); %& tpos ;
pos3=[find(pos1);find(pos2)];

Cgrm1=zeros(size(myCgrm));Cgrm1(:,pos1)=myCgrm(:,pos1);
Cgrm2=zeros(size(myCgrm));Cgrm2(:,pos2)=myCgrm(:,pos2);
%%

%%
figure(1);clf
ISDRAWCLUST=true;
NAMS={'LL','LH','RL','RH','BL','BH'};
NAMS={'Left-Low','Left-High','Right-Low','Right-High','Binaural-Low','Binaural-high'};


%subplot(2,1,1);
figure(4);clf;
%figure(2);clf
%imagesc((1:size(myCgrm,1))/fs,1:NFB3,log(abs((myCgrm'+eps))));title('spectrogram');axis xy;hold on;
imagesc((1:size(myCgrm,1))/fs,1:(NFB3-1),sqrt(abs((myCgrm(:,1:(NFB3-1))'+eps))));title('spectrogram');axis xy;hold on;
plot([0,size(myCgrm,1)],[NFB,NFB],'y-','LineWidth',3);
plot([0,size(myCgrm,1)],[2*NFB,2*NFB],'y-','LineWidth',3);
set(gca,'Ytick',[14,22,14+NFB,22+NFB,14+NFB*2,22+NFB*2]);
set(gca,'YtickLabel',NAMS);
xlabel('time (sec)');
if ISDRAWCLUST
    whr=[14, 22,14+NFB,22+NFB,14+NFB*2,22+NFB*2];
    whc=pos1(whr);
    CLRS={'b-','r-'};
    
    for I=1:length(whr)
        plot([0,size(myCgrm,1)],[whr(I)-2,whr(I)-2],CLRS{whc(I)+1},'LineWidth',2);
        plot([0,size(myCgrm,1)],[whr(I)+2,whr(I)+2],CLRS{whc(I)+1},'LineWidth',2);
        
    end
end
colormap(bone);
%subplot(2,1,2);
figure(3);clf

colormap(bone);
nm=(myhist+eps);


imagesc(nm);title('temporal coherence matrix');axis xy;hold on;

plot([1 NFB3],[NFB,NFB],'y-','LineWidth',3);
plot([NFB,NFB],[1 NFB3],'y-','LineWidth',3);
plot([1 NFB3],[2*NFB,2*NFB],'y-','LineWidth',3);
plot([2*NFB,2*NFB],[1 NFB3],'y-','LineWidth',3);
plot([1 NFB3],[NFB,NFB],'y-','LineWidth',3);
plot([NFB,NFB],[1 NFB3],'y-','LineWidth',3);
plot([1 NFB3],[2*NFB,2*NFB],'y-','LineWidth',3);
plot([2*NFB,2*NFB],[1 NFB3],'y-','LineWidth',3);
set(gca,'Xtick',[14,22,14+NFB,22+NFB,14+NFB*2,22+NFB*2]);
set(gca,'XtickLabel',NAMS);
set(gca,'Ytick',[14,22,14+NFB,22+NFB,14+NFB*2,22+NFB*2]);
set(gca,'YtickLabel',NAMS);

if ISDRAWCLUST
    whr=[14, 22,14+NFB,22+NFB,14+NFB*2,22+NFB*2];
    whc=pos1(whr);
    CLRS={'b-','r-'};
    
    for I=1:length(whr)
        plot([0,NFB3],[whr(I)-2,whr(I)-2],CLRS{whc(I)+1},'LineWidth',1);
        plot([0,NFB3],[whr(I)+2,whr(I)+2],CLRS{whc(I)+1},'LineWidth',1);
        plot([whr(I)-2,whr(I)-2],[0,NFB3],CLRS{whc(I)+1},'LineWidth',1);
        plot([whr(I)+2,whr(I)+2],[0,NFB3],CLRS{whc(I)+1},'LineWidth',1);
    end
end




%%

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


%     audiowrite('SPEC-mix.wav',0.1*ts,fs);
%     audiowrite('SPEC-unmix1.wav',0.1*Y1,fs);
%     audiowrite('SPEC-unmix2.wav',0.1*Y2,fs);
