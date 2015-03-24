clear all;close all;clc
addpath('/Users/jacoby/Dropbox (PPCA)/Research MIT/toolboxes/Sound_Texture_Synthesis_Toolbox');
[y,fs]=audioread('Wind.wav');
%[y,fs]=audioread('Insects.wav');
y=0.03*y/sum(y.^2);

nori_doplay(y,fs);
S0 = nori_generate_modsubands(y, fs);

y1=randn(size(y));
y1=0.03*y1/sum(y1.^2);
ITER=1;
for k=1:ITER,
    S1 = nori_generate_modsubands(y1, fs);
    
    new_mod_subbands=(S1.all_mod_subbands./S1.all_mod_subbands_envs).*S0.all_mod_subbands_envs;
    %new_mod_subbands=S0.all_mod_subbands;
    
    [collapse_subband_envs_n_up,collapse_subband_envs_n]=nori_collapse_modsubands_to_subbands(new_mod_subbands,S1.mod_filts,S1.P.audio_sr,S1.P.env_sr);

    collapse_subbands_1=(S1.subbands./(S1.subband_envs_up.^(1/S0.P.comp_exponent))).*(collapse_subband_envs_n_up.^(1/S0.P.comp_exponent));
   
    collapse_sound_1=collapse_subbands(collapse_subbands_1,S1.audio_filts);
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
    
    mysound=z;
    p2 = audioplayer(mysound/max(mysound), fs);p2.play;
   
end


mysound=z;
p2 = audioplayer(mysound/max(mysound), fs);p2.play

  