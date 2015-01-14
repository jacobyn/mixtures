clear all;close all;clc;
disp(sprintf('reading and formatting data...'));
%  load('FEATURES-Moz-v2.mat'); %high resolution
load('FEATURES-N-v4.mat');

targetfile1='C:\Users\user\Dropbox (PPCA)\Research MIT\mixtures\timit-test\si458.wav';
targetfile2='C:\Users\user\Dropbox (PPCA)\Research MIT\mixtures\timit-test\si664.wav';

T0=(300/1000*fs);T1=T0+200/1000*fs;

[ts1, fs]=audioread(targetfile1,[T0,T1]);
ts1=ts1/sqrt(mean(ts1.^2))*0.03; %normalize wav with rms
[ts2, fs]=audioread(targetfile2,[T0,T1]);
ts2=ts2/sqrt(mean(ts2.^2))*0.03; %normalize wav with rms
nn=min(length(ts1),length(ts2));
ts=ts1(1:nn)+ts2(1:nn);

KDT=KDTreeSearcher(FEATURES{1});

disp(sprintf('resynth...'));

[frames,sol_idxs,n3]=mozaic_make_frames(ts,NTP,fs,NFB,LF,HF,JMP);
sol=nan(length(sol_idxs),1);
for III=1:length(sol_idxs),
    II=sol_idxs(III);
    snapr=frames(III,:);
    
    [IDX,D] = knnsearch(KDT,snapr);
    sol(III,1)=IDX;
    
    mysnap=FEATURES{1}(IDX,:);
    fprintf('%d of %d %% %3.3g [ D= %g ]\n',III,length(sol_idxs),100*III/length(sol_idxs),D);
end


[frames3a,sol_idxs,n3]=mozaic_make_frames(ts1,NTP,fs,NFB,LF,HF,JMP);
sola=nan(length(sol_idxs),1);
for III=1:length(sol_idxs),
    II=sol_idxs(III);
    snapr=frames3a(III,:);
    
    [IDX,D] = knnsearch(KDT,snapr);
    sola(III,1)=IDX;
    
    mysnap=FEATURES{1}(IDX,:);
    fprintf('%d of %d %% %3.3g [ D= %g ]\n',III,length(sol_idxs),100*III/length(sol_idxs),D);
end


[frames3b,sol_idxs,n3]=mozaic_make_frames(ts2,NTP,fs,NFB,LF,HF,JMP);
solb=nan(length(sol_idxs),1);
for III=1:length(sol_idxs),
    II=sol_idxs(III);
    snapr=frames3b(III,:);
    
    [IDX,D] = knnsearch(KDT,snapr);
    solb(III,1)=IDX;
    
    mysnap=FEATURES{1}(IDX,:);
    fprintf('%d of %d %% %3.3g [ D= %g ]\n',III,length(sol_idxs),100*III/length(sol_idxs),D);
end
playme3a=mozaic_synth_sol(sola,INFO{1});p=audioplayer(ts1*10,fs);p.play;pause(1);p=audioplayer(playme3a*10,fs);p.play;pause(1)
playme3b=mozaic_synth_sol(solb,INFO{1});p=audioplayer(ts2*10,fs);p.play;pause(1);p=audioplayer(playme3b*10,fs);p.play;pause(1)

playme3=mozaic_synth_sol2([sola,solb],INFO{1});
[frames3,sol_idxs,n3]=mozaic_make_frames(playme3,NTP,fs,NFB,LF,HF,JMP);
nothing=zeros(.5*fs,1);
%%

sol2=[sol,sol];
oscore=inf;
%%
NITER=10000;
for KKK=1:NITER
    idx=randi(length(sol2));
    idxab=randi(2);
    idxw=randi(size(FEATURES{1},1));
    nsol2=sol2;
    nsol2(idx,idxab)=idxw;
    playme2=mozaic_synth_sol2(nsol2,INFO{1});
    % p=audioplayer(playme2*10,fs);p.play;
    [frames2,sol_idxs,n3]=mozaic_make_frames(playme2,NTP,fs,NFB,LF,HF,JMP);
    
    playme2a=mozaic_synth_sol(nsol2(:,1),INFO{1});
    playme2b=mozaic_synth_sol(nsol2(:,2),INFO{1});
    [frames2a,~,~]=mozaic_make_frames(playme2a,NTP,fs,NFB,LF,HF,JMP);
    [frames2b,~,~]=mozaic_make_frames(playme2b,NTP,fs,NFB,LF,HF,JMP);
    
    
    
    NFF=min(size(frames,1),size(frames2,1));
    % figure(1);clf;imagesc(corr(frames(NFF,:)',frames2(NFF,:)'));
    score2=(sqrt(sum(sum( (frames(1:NFF,:)-frames2(1:NFF,:)).^2))))/1000;
    score2ab=sum(abs(diag(corr(frames2a',frames2b'))'));
     score=score2+score2ab/10;


    score3=(sqrt(sum(sum( (frames(1:NFF,:)-frames3(1:NFF,:)).^2))))/1000;
    score3ab=sum(abs(diag(corr(frames3a',frames3b'))'));
    score3t=score3+score3ab/10;

     
    if score<=oscore || (mod(KKK,100)==1)
        sol2=nsol2;
        oscore=score;
         msg='';
%         for kk=1:2
%             for ll=1:size(sol2,1)
%                 msg=sprintf('%s %3d',msg,sol2(ll,kk));
%             end
%             msg=sprintf('%s\n',msg);
%         end
        fprintf('iter %d \t %g %g \t %g*** %g %g \t%g \n%s',KKK,score2,score2ab,score,score3,score3ab,score3t,msg);
        
        figure(1);clf;
        
        nf=1;
        for kk=1:nf,
            subplot(4,nf,kk);
            imagesc(mozaic_reorder_frames(frames3a,NFB,NTP));axis off;
            subplot(4,nf,kk+nf);
            imagesc(mozaic_reorder_frames(frames3b,NFB,NTP));axis off;
            subplot(4,nf,kk+nf*2);
            imagesc(mozaic_reorder_frames(frames2a,NFB,NTP));axis off;
            subplot(4,nf,kk+nf*3);
            imagesc(mozaic_reorder_frames(frames2b,NFB,NTP));axis off;
        end
        
         theplay=10*[playme2a;nothing;playme2b;nothing;nothing;playme3a;nothing;playme3b];
%         theplay=10*[playme2a;nothing;playme2b;nothing];
         p=audioplayer(theplay,fs);p.play;
%          p=audioplayer(playme3a*10,fs);p.play;pause(1);
%                   p=audioplayer(playme3b*10,fs);p.play;pause(1);
%                   pause(1);
%         p=audioplayer(playme2a*10,fs);p.play;pause(1);
% 
%         p=audioplayer(playme2b*10,fs);p.play;pause(1);
        
    else
        fprintf('iter %d \t %g %g\n',KKK,score2,score2ab);
    end
    
    
    
end






%         gsnap=max(snap.^2-(10.^(mysnap/10)),0);
%         grnap=10*log10((gsnap.*gsnap)+eps);
%         figure(1);clf;subplot(2,2,1);imagesc(reshape(snapr,NFB,NTP))
%         subplot(2,2,2);imagesc(reshape(mysnap,NFB,NTP))
%         subplot(2,2,3);imagesc(reshape(grnap,NFB,NTP))
