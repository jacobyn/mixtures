%clear all;close all;clc
cd '/Users/jacoby/Dropbox (PPCA)/Research MIT/mixtures'
load('WFR-200-MIX-WHATTYPE-A2.mat');
fs=FS;
S = nori_generate_fqsubands(soundM, fs);
audio_filts=S.audio_filts;

subband_envs=S.subband_envs;
measurement_win = ones(length(subband_envs),1);
mod_filt_length = length(subband_envs);%%
env_sr=S.P.env_sr;
N_mod_channels = 20; %These next four parameters control the modulation filterbank from which modulation power is measured
% P.low_mod_f = 0.5; %Hzl
low_mod_f = 1; %Hz %%%%%%%%%%NOTE NOTE NOTE!!! I HAD TO CHANGE IT BECAUSE I WORK WITH SHORT SAMPLES
hi_mod_f = 200; %Hz
mod_filt_Q_value = 2;
use_zp = 0;% 0 means circular convolution; 1 means zeropadding (for modulation filtering)

[mod_filts,Hz_mod_cfreqs,mod_freqs] = make_constQ_cos_filters(mod_filt_length, env_sr, N_mod_channels, low_mod_f, hi_mod_f, mod_filt_Q_value);%%



mod_power=nan(size(audio_filts,2),size(mod_filts,2));
for j=1:size(audio_filts,2) %go through subbands
    mod_power(j,:) = stat_mod_power_win(subband_envs(:,j), mod_filts, use_zp, measurement_win);
end

S = nori_measure_texture_stats_simplify(soundM, fs);

%% 
mymask=rand(size(subband_envs))>0.5;
subband_envs1=subband_envs.*mymask;
subband_envs2=subband_envs.*(~mymask);
mod_power1=nan(size(audio_filts,2),size(mod_filts,2));
mod_power2=nan(size(audio_filts,2),size(mod_filts,2));
for j=1:size(audio_filts,2) %go through subband
    mod_power1(j,:) = stat_mod_power_win(subband_envs1(:,j), mod_filts, use_zp, measurement_win);
    mod_power2(j,:) = stat_mod_power_win(subband_envs2(:,j), mod_filts, use_zp, measurement_win);
end

score=sum(sum(mod_power1.*wfr+mod_power2.*wfr));
optim=optimset('Display','iter');

Na=size(audio_filts,2);
Nb=size(mod_filts,2);
scores=nan(Na,1);
for j=1:size(audio_filts,2) %go through subband
    subband_envs_j=subband_envs(:,j);
    mask_j=mymask(:,j);
    wfr_j=wfr(j,:);
    
    X0=double(mask_j*1.0);
    mask_j_new = fmincon(@(x)temp_modband_line_score(x,wfr_j,subband_envs_j,Na,Nb,mod_filts,use_zp,measurement_win),X0,[],[],[],[],zerXos(size(X0)),ones(size(X0)),[],optim);
    scores(j)=temp_modband_line_score(mask_j_new,wfr_j,subband_envs_j,Na,Nb,mod_filts,use_zp,measurement_win);
end

