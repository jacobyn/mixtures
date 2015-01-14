%%
[ts,fs]=wavread('data/mix-v1.M.1.32.wav');

% generate filters
%[FilterBank,SubbandFrequencies]=make_erb_cos_filters(length(ts),fs,No_of_Subbands,Min_Freq,Max_Freq);

[FilterBank,SubbandFrequencies]=make_erb_cos_filters(length(ts),fs,50,100,3000);

% filter your time series into subbands
Cgrm= generate_subbands(ts,FilterBank);
%If you want to generate a time series from a Cochleagram use

figure(2);imagesc((1:length(ts))/fs,SubbandFrequencies,abs(Cgrm)');axis xy;ylabel('fq (Hz)');xlabel('time (s)');
%%

%%%NewTimeSeries=collapse_subbands(Cgrm,FilterBank)
[heart_scale_label, heart_scale_inst] = libsvmread('../heart_scale');
model = svmtrain(heart_scale_label, heart_scale_inst, '-c 1 -g 0.07');


[predict_label, accuracy, dec_values] = svmpredict(heart_scale_label, heart_scale_inst, model); % test the training data

% % get w and b
% w = model.SVs' * model.sv_coef;
% b = -model.rho;
% 
% if model.Label(1) == -1
%   w = -w;
%   b = -b;
% end

%%

logcs=-1:3;
loggs=-4:1;
%logcs=-5:.5:4;
%loggs=-10:.5:8;

mygrid=nan(length(logcs),length(loggs));


bestcv = 0;
for I=1:length(logcs),
  for J=1:length(loggs),
    log2c=logcs(I); 
    log2g=loggs(J);

    cmd = ['-v 5 -c ', num2str(2^log2c), ' -g ', num2str(2^log2g)];
    cv = svmtrain(heart_scale_label, heart_scale_inst, cmd);
    if (cv >= bestcv),
      bestcv = cv; bestc = 2^log2c; bestg = 2^log2g;
    end
    mygrid(I,J)=cv;
    fprintf('%g %g %g (best c=%g, g=%g, rate=%g)\n', log2c, log2g, cv, bestc, bestg, bestcv);
  end
end
figure(5);imagesc(logcs,loggs,mygrid);axis xy;