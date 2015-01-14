function [frames,sol_idxs,n3]=mozaic_make_frames(ts,NTP,fs,NFB,LF,HF,JMP)

[FilterBank,SubbandFrequencies]=make_erb_cos_filters(length(ts),fs,NFB-2,LF,HF); %make filter bank
Cgrm= generate_subbands(ts,FilterBank)';
Cgrm_jmp=Cgrm(:,1:JMP:end);
n3=round(NTP/3);
sol_idxs=1:n3:(size(Cgrm_jmp,2));

frames=[];
for III=1:length(sol_idxs),
    II=sol_idxs(III);
    T2=min(II+NTP-1,size(Cgrm_jmp,2));
    if (II+NTP-1)>size(Cgrm_jmp,2)
        snap=snap;
        snap(:,(1):(T2-II+1))=Cgrm_jmp(:,(II):(T2));  
    else
        snap=Cgrm_jmp(:,(II):(II+NTP-1));    
    end
    snapr=reshape(snap,1,NTP*NFB);
    snapr=10*log10((snapr.*snapr)+eps);
    if isempty(frames)
        frames=nan(length(sol_idxs),length(snapr));
        frames(1,:)=snapr;
    else
        frames(III,:)=snapr;
    end
end
         
