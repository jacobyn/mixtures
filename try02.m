[FilterBank,SubbandFrequencies]=make_erb_cos_filters(length(ts),fs,50,100,3000);
MM=1;
for KK=1:10,
    mfname=sprintf('data/mix-v1.M.%d.%d.wav',M,KK);

[ts,fs]=audioread(mfname);

% generate filters
%[FilterBank,SubbandFrequencies]=make_erb_cos_filters(length(ts),fs,No_of_Subbands,Min_Freq,Max_Freq);



% filter your time series into subbands
Cgrm= generate_subbands(ts,FilterBank);
%If you want to generate a time series from a Cochleagram use


figure(1);
subplot(2,2,1);
imagesc((1:length(ts))/fs,SubbandFrequencies,abs(Cgrm)');axis xy;ylabel('fq (Hz)');xlabel('time (s)');


NFB=length(SubbandFrequencies);
myhist=zeros(NFB,NFB);

myCgrm=abs(Cgrm)/max(max(abs(Cgrm)));
for I=1:NFB,
    for J=1:NFB,
         if I~=J
        myhist(I,J)=sum( (myCgrm(:,I)).*( myCgrm(:,J)));
         else
             myhist(I,J)=sum( (myCgrm(:,I)).*( myCgrm(:,J)));
         end
    end
end
subplot(2,2,2);
imagesc(SubbandFrequencies,SubbandFrequencies,myhist);

 AP=audioplayer(ts,fs);
AP.play