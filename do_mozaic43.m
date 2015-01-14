clear all;close all;clc;
disp(sprintf('reading and formatting data...'));
%  load('FEATURES-Moz-v2.mat'); %high resolution
 load('FEATURES-N-v4.mat');

targetfile='C:\Users\user\Dropbox (PPCA)\Research MIT\mixtures\timit-test\si458.wav';
targetfile='C:\Users\user\Dropbox (PPCA)\Research MIT\mixtures\timit-test\si664.wav';
[ts, fs]=audioread(targetfile);
ts=ts/sqrt(mean(ts.^2))*0.03; %normalize wav with rms
[FilterBank,SubbandFrequencies]=make_erb_cos_filters(length(ts),fs,NFB-2,LF,HF); %make filter bank
Cgrm= generate_subbands(ts,FilterBank)';
Cgrm_jmp=Cgrm(:,1:JMP:end);
n3=round(NTP/3);
KDT=KDTreeSearcher(FEATURES{1});

info=cell(size(FEATURES{1},1),1);cnt=0;
for J=1:size(FEATURES{1},1),
    cnt=cnt+1;
    info{cnt,1}.fname=INFO{1}{J}.fname;
    info{cnt,1}.range=INFO{1}{J}.range;
end
disp(sprintf('resynth...'));
playme=zeros(1000,1);
for II=1:n3:(size(Cgrm_jmp,2)-NTP-1),
    
    snap=Cgrm_jmp(:,(II):(II+NTP-1));
    snap=reshape(snap,1,NTP*NFB);
    snapr=10*log10((snap.*snap)+eps);
    
%     for JJ=1:size(FEATURES{1},1)
%         asnap=FEATURES{1}(J,:);
         
    [IDX,D] = knnsearch(KDT,snapr);
   mysnap=FEATURES{1}(IDX,:);
        gsnap=max(snap.^2-(10.^(mysnap/10)),0);
        grnap=10*log10((gsnap.*gsnap)+eps);
%         figure(1);clf;subplot(2,2,1);imagesc(reshape(snapr,NFB,NTP))
%         subplot(2,2,2);imagesc(reshape(mysnap,NFB,NTP))
%         subplot(2,2,3);imagesc(reshape(grnap,NFB,NTP))

    
    fprintf('%d of %d %% %3.3g [ D= %g ]\n',II,size(Cgrm_jmp,2),100*II/size(Cgrm_jmp,2),D);
    [Y, fs]=audioread(info{IDX}.fname, info{IDX}.range);
    Yh=Y.*hanning(length(Y));
    n3=round(length(Yh)/3);
    playmeold=playme;
    playme=zeros(length(playmeold)+n3,1);
    playme(1:length(playmeold))=playmeold;
    playme((length(playmeold)-(2*n3)+1):(length(playmeold)-(2*n3)+length(Yh)))=playme((length(playmeold)-(2*n3)+1):(length(playmeold)-(2*n3)+length(Yh)))+Yh;
    
end
%%
p=audioplayer(10*ts,fs);p.play;
%%
p=audioplayer(playme*10,fs);p.play;
%%
% % % % % n2t=NTP/2
% % % % % playme=zeros(1000,1);
% % % % % for II=1:n2t:(size(Cgrm_jmp,2)-NTP-1),
% % % % %     
% % % % %     snap=Cgrm_jmp(:,(II):(II+NTP-1));
% % % % %     snapr=reshape(snap,1,NTP*NFB);
% % % % %     snapr=10*log10((snapr.*snapr)+eps);
% % % % % 
% % % % %     [IDX,D] = knnsearch(KDT,snapr);
% % % % %     fprintf('%d of %d %% %3.3g [ D= %g ]\n',II,size(Cgrm_jmp,2),100*II/size(Cgrm_jmp,2),D);
% % % % %     [Y, fs]=audioread(info{IDX}.fname, info{IDX}.range);
% % % % %     Yh=Y.*hanning(length(Y));
% % % % %     n2=round(length(Yh)/2);
% % % % %     playmeold=playme;
% % % % %     playme=zeros(length(playmeold)+n2,1);
% % % % %     playme(1:length(playmeold))=playmeold;
% % % % %     playme((length(playmeold)-(n2)+1):(length(playmeold)-(n2)+length(Yh)))=playme((length(playmeold)-(n2)+1):(length(playmeold)-(n2)+length(Yh)))+Yh;
% % % % %     
% % % % % end
% % % % % %%
% % % % % p=audioplayer(playme,fs);p.play;
