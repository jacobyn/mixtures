clear all;close all;clc
addpath('/Users/jacoby/Dropbox (PPCA)/Research MIT/toolboxes/Sound_Texture_Synthesis_Toolbox');
[y,fs]=audioread('Wind.wav');
y=0.03*y/sum(y.^2);

%nori_doplay(y,fs);
S0 = nori_generate_modsubands(y, fs);

%[collapse_subband_envs_0_up,collapse_subband_envs_0]=nori_collapse_modsubands_to_subbands(S0);
%S.subband_envs,S.all_mod_subbands,S.P.audio_sr,S.P.env_sr;size(S.subband_envs)
[collapse_subband_envs_0_up,collapse_subband_envs_0]=nori_collapse_modsubands_to_subbands(S0.all_mod_subbands,S0.mod_filts,S0.P.audio_sr,S0.P.env_sr);


% %collapse mod-subbands to subbands
% collapse_subband_envs_0=zeros(size(S0.subband_envs));
% for j=1:size(S0.audio_filts,2) %go through subbands
%     my_mod_subbands=reshape(S0.all_mod_subbands(j,:,:),size(S0.all_mod_subbands,2),size(S0.all_mod_subbands,3));
%     my_colapse_envs=collapse_subbands(my_mod_subbands, S0.mod_filts);
%     my_colapse_envs(my_colapse_envs<0)=0;
%     collapse_subband_envs_0(:,j)=my_colapse_envs;
% end
% collapse_subband_envs_0_up=resample(collapse_subband_envs_0,S0.P.audio_sr,S0.P.env_sr);

y1=randn(size(y));
y1=0.03*y1/sum(y1.^2);
ITER=5;
for k=1:ITER,
    S1 = nori_generate_modsubands(y1, fs);
    orig_subband_envs_1_up=resample(S1.subband_envs,S1.P.audio_sr,S1.P.env_sr);
    collapse_subbands_1=(S1.subbands./(orig_subband_envs_1_up+eps)).*collapse_subband_envs_0_up;
    collapse_sound_1=collapse_subbands(collapse_subbands_1,S0.audio_filts);
    z=collapse_sound_1;
    y1=z;
    
    
    
    
    S3 = nori_generate_modsubands(z, fs);
    S0s = nori_measure_texture_stats_simplify(y, fs);
    S3s = nori_measure_texture_stats_simplify(z, fs);
    
    K=3;J=12;
    figure(1);clf;
    plot((S0.all_mod_subbands_envs(K,:,J)));hold all;
    plot(S1.all_mod_subbands_envs(K,:,J));hold all;
    plot(S3.all_mod_subbands_envs(K,:,J));hold all;
    
    figure(2);
    subplot(2,2,1);
    imagesc(S0s.mod_power');axis xy;
    
    subplot(2,2,2);
    imagesc(S3s.mod_power');axis xy;
    
    subplot(2,2,3);
    imagesc(S0s.subband_envs');axis xy;
    
    subplot(2,2,4);
    imagesc(S3s.subband_envs');axis xy;
    drawnow;
    figure(3);clf;
    plot(S0s.env_mean);hold all;
    plot(S3s.env_mean);
    
    nori_doplay(z,fs);
end




