function rframes=mozaic_reorder_frames(frames,NFB, NTP)
rframes=[];
for kk=1:size(frames,1)
    mframe=reshape(frames(kk,:),[NFB NTP]);
    rframes=[rframes,mframe];    
end
    
end