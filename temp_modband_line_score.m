function score=temp_modband_line_score(mask_j,wfr_j,subband_envs_j,Na,Nb,mod_filts,use_zp,measurement_win)
% Na=size(audio_filts,2);
% Nb=size(mod_filts,2);
mod_power1=nan(Na,Nb);
mod_power2=nan(Na,Nb);
subband_envs_j1=subband_envs_j.*(mask_j);
subband_envs_j2=subband_envs_j.*(1-mask_j);

for j=1:Na %go through subband
    mod_power1 = stat_mod_power_win(subband_envs_j1, mod_filts, use_zp, measurement_win);
    mod_power2 = stat_mod_power_win(subband_envs_j2, mod_filts, use_zp, measurement_win);
end

score=sum(sum(mod_power1.*wfr_j+mod_power2.*wfr_j));
end