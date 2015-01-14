 close all;clc;clear all;


NFB=50;
Ms=9:-1:1;
ITER=10;

for M=Ms,
    myhist_all=zeros(NFB,NFB);
    
    
    for KK=1:ITER,
        mfname=sprintf('data/mix-v1.M.%d.%d.wav',M,KK);
        display(mfname);
        
        [ts,fs]=audioread(mfname);
        [FilterBank,SubbandFrequencies]=make_erb_cos_filters(length(ts),fs,NFB-2,100,3000);
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
        nyhist=myhist/max(max(myhist));
        myhist_all=myhist_all+myhist;
        
        figure(1);
        subplot(2,2,1);
        imagesc((1:length(ts))/fs,SubbandFrequencies,abs(Cgrm)');axis xy;ylabel('fq (Hz)');xlabel('time (s)');
        colorbar;
        
        figure(1);
        subplot(2,2,2);
        imagesc(SubbandFrequencies,SubbandFrequencies,myhist);
        colorbar;
        
        subplot(2,2,3);
        imagesc(SubbandFrequencies,SubbandFrequencies,myhist_all);
        colorbar;
        drawnow
        pause;
        %     AP=audioplayer(ts,fs);
        %     AP.play
        %
        
    end
    figure(2);
    subplot(3,3,M);
    imagesc(SubbandFrequencies,SubbandFrequencies,myhist_all);
    title(sprintf('num spks %d',M));
end
%%



%%
% logcs=-5:.5:4;
% loggs=-10:.5:8;
% 
% mygrid=nan(length(logcs),length(loggs));
% 
% 
% bestcv = 0;
% for I=1:length(logcs),
%   for J=1:length(loggs),
%     log2c=logcs(I); 
%     log2g=loggs(J);
% 
%     cmd = ['-v 5 -c ', num2str(2^log2c), ' -g ', num2str(2^log2g)];
%     cv = svmtrain(labels, vals, cmd);
%     if (cv >= bestcv),
%       bestcv = cv; bestc = 2^log2c; bestg = 2^log2g;
%     end
%     mygrid(I,J)=cv;
%     fprintf('%g %g %g (best c=%g, g=%g, rate=%g)\n', log2c, log2g, cv, bestc, bestg, bestcv);
%   end
% end
