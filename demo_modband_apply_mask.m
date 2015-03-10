% this is doing the masking on mod sub bands!
S0 = nori_generate_modsubands(soundM, FS0);
S1 = nori_generate_modsubands(mixme{1}, FS0);
S2 = nori_generate_modsubands(mixme{2}, FS0);

%nmask=S1.all_mod_subbands_envs>S2.all_mod_subbands_envs;
mmask=S1.all_mod_subbands_envs_dn>S2.all_mod_subbands_envs_dn;
nmask=zeros(size(S0.all_mod_subbands));
for J=1:size(mmask,1)
    for I=1:size(mmask,3)
        myw=double(mmask(J,:,I));
        nmask(J,:,I)=resample(myw(:),S0.P.env_sr,S0.P.mod_env_sr);
    end
end
%%
%mask=resample(double(mask),S0.P.env_sr,S0.P.mod_env_sr);

%mask=S1.all_mod_subbands_envs>S2.all_mod_subbands_envs;

%mask(:,:,1)=0.5;

 collapse_subband_envs_1=zeros(size(S0.subband_envs));
for j=1:size(S0.audio_filts,2) %go through subbands
     my_mod_subbands=reshape(S0.all_mod_subbands(j,:,:).*nmask(j,:,:),size(S0.all_mod_subbands,2),size(S0.all_mod_subbands,3));
     my_colapse_envs=collapse_subbands(my_mod_subbands, S0.mod_filts);
     %my_colapse_envs=abs(my_colapse_envs);
     my_colapse_envs(my_colapse_envs<0)=0;
     %my_colapse_envs=my_colapse_envs.^(S0.P.comp_exponent);
     collapse_subband_envs_1(:,j)=my_colapse_envs;
end

 collapse_subband_envs_2=zeros(size(S0.subband_envs));
for j=1:size(S0.audio_filts,2) %go through subbands
     my_mod_subbands=reshape(S0.all_mod_subbands(j,:,:).*(1-nmask(j,:,:)),size(S0.all_mod_subbands,2),size(S0.all_mod_subbands,3));
     my_colapse_envs=collapse_subbands(my_mod_subbands, S0.mod_filts);
     %my_colapse_envs=abs(my_colapse_envs);
     my_colapse_envs(my_colapse_envs<0)=0;
     %my_colapse_envs=my_colapse_envs.^(S0.P.comp_exponent);
     collapse_subband_envs_2(:,j)=my_colapse_envs;
end



collapse_subband_envs_1_up=resample(collapse_subband_envs_1,S0.P.audio_sr,S0.P.env_sr);

%subband_envs = subband_envs.^P.comp_exponent;
collapse_subbands_1=S0.subbands.*collapse_subband_envs_1_up;
collapse_sound_1=collapse_subbands(collapse_subbands_1,S0.audio_filts);

collapse_subband_envs_2_up=resample(collapse_subband_envs_2,S0.P.audio_sr,S0.P.env_sr);
collapse_subbands_2=S0.subbands.*collapse_subband_envs_2_up;
collapse_sound_2=collapse_subbands(collapse_subbands_2,S0.audio_filts);


     
        p1 = audioplayer(collapse_sound_1/max(collapse_sound_1), FS);p1.play
        pause(1.1*length(collapse_sound_1)/FS +0.5);
        
        p2 = audioplayer(collapse_sound_2/max(collapse_sound_2), FS);p2.play
        pause(1.1*length(collapse_sound_2)/FS +0.5);
        
        p2 = audioplayer(soundM/max(soundM), FS);p2.play
        pause(1.1*length(soundM)/FS +0.5);
    %%    
%  
%  figure(11);clf;    
%  for j=1:size(S1.audio_filts,2) %go through subbands
%      subplot(6,6,j);
%      imagesc(log2(abs(reshape(S1.all_mod_subbands(j,2:end,2:end),size(S1.all_mod_subbands,2)-1,size(S1.all_mod_subbands,3)-1)))')
%  end
%  
%  figure(12);clf;    
%  for j=1:size(S1.audio_filts,2) %go through subbands
%      subplot(6,6,j);
%      imagesc(log2(abs(reshape(S2.all_mod_subbands(j,2:end,2:end),size(S2.all_mod_subbands,2)-1,size(S2.all_mod_subbands,3)-1)))')
%  end
%    
  %%
 figure(13);clf;    
 for j=1:min(size(S1.audio_filts,2),4) %go through subbands
     subplot(2,2,j);
     plot(S0.subband_envs(:,end-j),'k-');hold on;
    
     plot(S2.subband_envs(:,end-j),'r-');hold on;
      plot(S1.subband_envs(:,end-j),'b-');hold on;
     plot(collapse_subband_envs_1(:,(end-j)),'c+-','LineWidth',2.5);hold on;
 end
        
%%
% this is doing the masking on sub bands (directly without using mod-subands info)!

S0 = nori_generate_modsubands(soundM, FS0);
S1 = nori_generate_modsubands(mixme{1}, FS0);
S2 = nori_generate_modsubands(mixme{2}, FS0);
mask=S1.subband_envs>S2.subband_envs;


collapse_subband_envs_1=S0.subband_envs.*mask;
collapse_subband_envs_2=S0.subband_envs.*(~mask);

collapse_subband_envs_1_up=resample(collapse_subband_envs_1,S0.P.audio_sr,S0.P.env_sr);
collapse_subbands_1=S0.subbands.*collapse_subband_envs_1_up;
collapse_sound_1=collapse_subbands(collapse_subbands_1,S0.audio_filts);

collapse_subband_envs_2_up=resample(collapse_subband_envs_2,S0.P.audio_sr,S0.P.env_sr);
collapse_subbands_2=S0.subbands.*collapse_subband_envs_2_up;
collapse_sound_2=collapse_subbands(collapse_subbands_2,S0.audio_filts);


     
        p1 = audioplayer(collapse_sound_1/max(collapse_sound_1), FS);p1.play
        pause(1.1*length(collapse_sound_1)/FS +0.5);
        
        p2 = audioplayer(collapse_sound_2/max(collapse_sound_2), FS);p2.play
        pause(1.1*length(collapse_sound_2)/FS +0.5);
        
        %p2 = audioplayer(soundM/max(soundM), FS);p2.play
        %pause(1.1*length(soundM)/FS +0.5);
  
        %%
  figure(13);      
  for j=1:min(size(S1.audio_filts,2),4) %go through subbands
    subplot(2,2,j);
     plot(collapse_subband_envs_1(:,j),'g:','LineWidth',1.5);hold on;
 end
             
