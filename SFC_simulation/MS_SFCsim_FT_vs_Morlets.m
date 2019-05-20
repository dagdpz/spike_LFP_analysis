% Compare the performance of Fourier transform and Complex Morlet Wavelet
% convolution

clc, clear, close all;
% 
%     keys.path_to_save = sprintf('%s\\ft_vs_morlet', pwd);
% 
%     if exist(sprintf('%s\\ft_vs_morlet', pwd)) == 0
%         mkdir(sprintf('%s\\ft_vs_morlet', pwd))
%     end
% 
%     FigName_short1= sprintf('FT_vs_Morlet'); %name of the file
%     FigName_Long1=sprintf('FT convolution vs Morlet wavelet convolution'); % overall title of the figure
%     set(0,'DefaultTextInterpreter','none');
%     wanted_papersize= [21 6];
%     h1=figure('outerposition',[0 0 1900 1200],'name', FigName_Long1);
%     figure(h1)
num_trials = [5, 10, 50, 100, 500, 1000]; %
%time_cmp = zeros(numel(num_trials), 4);
for numtr = 1:numel(num_trials);
    
    iter = 1;
    while iter <= 30,
    
    
        cfg = [];
        cfg.numtrl	= num_trials(numtr);
        cfg.fsample     = 1000; % Hz
        cfg.trllen      = 5; % s

        cfg.s1.freq     = 10; % Hz
        cfg.s1.phase    = 'random'; % 'random' or 0
        cfg.s1.ampl     = 1;

        cfg.s2.freq = 25; % Hz
        cfg.s2.ampl = 0.5;
        cfg.s2.phase = 0;

        cfg.s3.freq = 60; % Hz
        cfg.s3.ampl = 0.2;
        cfg.s3.phase = 0;

        cfg.noise.ampl = 0.2;

        cfg.method	= 'superimposed';
        cfg.output	= 'all'; % mixed or all

        data = ft_freqsimulation(cfg);

        lfp_freq = [cfg.s1.freq cfg.s2.freq cfg.s3.freq];
        lfp_amp =  [cfg.s1.ampl cfg.s2.ampl cfg.s3.ampl];

        % spikes
        spikeRate = [40]; % Hz, overall spike rate
        spikePhaseFreq = 1; % 1 or 2 or 3, s1 or s2 or s3 components
        spikePhaseMean = [-pi]; % align to trough
        spikePhaseStd  = [pi/60]; % rad, spread around spikePhaseMean
        spikeLockProb = 0.5; % probability of spike to be locked, the phase-locked firing rate would be max(spikeRate*spikeLockProb,lfp_freq(spikePhaseFreq))

        lockedSpikes = zeros(cfg.numtrl,cfg.trllen*cfg.fsample);
        Spikes = MS_test_FT_simulate_spike_train(spikeRate*(1-spikeLockProb),cfg.trllen,cfg.numtrl); % non-locked spikes
        for t = 1:cfg.numtrl,
            % for each trial, find the phase (of the trough) of the corresponding freq component s
             s = data.trial{t}(spikePhaseFreq+1,:);
             taxis = data.time{t};
             p = angle(hilbert(s));
             trough_idx = find(diff(p)<(min(p)+abs(min(p))/20)); % find indices of the troughs

             lockedSpikes_idx = trough_idx(rand(size(trough_idx)) <= spikeRate*spikeLockProb/lfp_freq(spikePhaseFreq));

             % add some variability to relative phase of spikes and troughs
        % 	 lockedSpikes_idx = fix(lockedSpikes_idx+randn(size(lockedSpikes_idx))*spikePhaseStd*(cfg.fsample/lfp_freq(spikePhaseFreq)));
        % 	 lockedSpikes_idx = ig_limit_range_min_max(lockedSpikes_idx,1,cfg.fsample*cfg.trllen);
             lockedSpikes(t,lockedSpikes_idx)=1;
             Spikes(t,lockedSpikes_idx) = 1; % add locked spikes to spikes
             data.trial{t} = [data.trial{t}; Spikes(t,:)]; 

             if 0 % plot lfp and spikes (for illustration and debugging)
                 spikes_idx = find(Spikes(t,:));
                 rand_rgb = rand(1,3);
                 plot(taxis,s,'Color',rand_rgb); hold on;
                 plot(taxis(spikes_idx),-1*lfp_amp(spikePhaseFreq),'rx','MarkerEdgeColor',rand_rgb); % all spikes
                 if ~isempty(lockedSpikes_idx)
                    plot(taxis(lockedSpikes_idx),-1*lfp_amp(spikePhaseFreq),'ro','MarkerEdgeColor',rand_rgb); % locked spikes
                 end
                 pause;
             end

        end

        disp(sprintf('spikeRate %d Hz, spikeLockProb %.2f, mean rate %.2f, mean locked rate %.2f %d spikes',spikeRate,spikeLockProb,mean(sum(Spikes,2)/cfg.trllen),mean(sum(lockedSpikes,2)/cfg.trllen),sum(sum(Spikes))));
        data.label = {'lfp1', 's1', 's2', 's3', 'noise', 'spk1'};

        %% MTMCONVOL
        disp('Spike-triggered spectrum, mtmconvol');
        tic
        cfg           = [];
        cfg.method    = 'mtmconvol';
        cfg.foi       = 10.^(0.6:0.05:2);
        cfg.taper     = 'hanning';
        cfg.spikechannel = data.label{6};
        cfg.channel      = data.label{1};

        cfg.borderspikes='no'; % LS 
        cfg.t_ftimwin = 6./cfg.foi; % X cycles per frequency %repmat(0.5,size(cfg.foi));% fixed windows, LS 

        stsConvol     = ft_spiketriggeredspectrum(cfg, data);
        ft_time = toc;

    %     % compute the statistics on the phases
    %     cfg               = [];
    %     cfg.method        = 'ppc1'; % compute the Pairwise Phase Consistency, can also be plv, ang, ppc1, ppc2, ral
    %     cfg.spikechannel  = stsConvol.label{1};
    %     cfg.channel       = stsConvol.lfplabel(1); % selected LFP channels
    %     cfg.avgoverchan   = 'unweighted'; % weight spike-LFP phases irrespective of LFP power
    %     cfg.timwin        = 'all'; % compute over all available spikes in the window
    %     cfg.latency       = 'maxperiod'; % [0 nanmax(stsConvol.trialtime(:))]; %
    %     statSts           = ft_spiketriggeredspectrum_stat(cfg,stsConvol);

        %% MORLETS

        keys.spike_field               = [];
        keys.spike_field.method        = 'ppc1'; % compute the Pairwise Phase Consistency, can also be plv, ang, ppc1, ppc2, ral
        keys.spike_field.spikechannel  = 'spk1';
        keys.spike_field.channel       = 'lfp1'; % selected LFP channels
        keys.spike_field.avgoverchan   = 'unweighted'; % weight spike-LFP phases irrespective of LFP power
        keys.spike_field.timwin        = 'all'; % compute over all available spikes in the window
        keys.spike_field.latency       = 'maxperiod'; % [0 nanmax(stsConvol.trialtime(:))]; %
        keys.spike_field.N_cycles      = 6;

        tic
        keys.LFP.frequencies        =10.^(0.6:0.05:2);
        keys.LFP.binsize            =1.1*keys.spike_field.N_cycles./keys.LFP.frequencies;
        keys.LFP.Morlet_bandwidth   =2*(keys.spike_field.N_cycles./(2*pi*keys.LFP.frequencies)).^2;
        keys.LFP.SR                 =1000; %1017.252604166667;

        %morlet_borders=ceil(keys.LFP.binsize/2*keys.LFP.SR)/keys.LFP.SR;
        for f=1:numel(keys.LFP.frequencies)
            fr=keys.LFP.frequencies(f);
            morlet_borders=ceil(keys.LFP.binsize(f)/2*keys.LFP.SR)/keys.LFP.SR;
            keys.LFP.Morlet_Matrix{f}=cmorwavf(-morlet_borders,morlet_borders,2*morlet_borders*keys.LFP.SR,keys.LFP.Morlet_bandwidth(f),fr);
        end

        % Morlet wavelets convlution
        nspikesaccum=0;

        for t=1:numel(data.trial)
            clear spectrum
            for f=1:numel(keys.LFP.Morlet_Matrix)
                filtered(f,:)=conv(data.trial{t}(1,:),keys.LFP.Morlet_Matrix{f},'same');
            end
            spectrum(:,:,1)=filtered(:,data.trial{t}(end,:)==1);
            nspikesend= nspikesaccum+size(spectrum,2);
            spectrum_perm(nspikesaccum+1:nspikesend,1,:)=permute(spectrum,[2,3,1]);
            timevec(nspikesaccum+1:nspikesend,1) = find(data.trial{t}(end,:)==1)/1000;
            trialvec(nspikesaccum+1:nspikesend,1) = t;
            trialtimevec(t, :) = [min(data.time{t}), max(data.time{t})] ;
            nspikesaccum=nspikesend;
        end
        morlet_time = toc;

    %     morlet.lfplabel = {keys.spike_field.channel};
    %     morlet.label = {keys.spike_field.spikechannel};
    %     morlet.freq = keys.LFP.frequencies;
    %     morlet.fourierspctrm = {spectrum_perm}; 
    %     morlet.time = {timevec};
    %     morlet.trial = {trialvec};
    %     morlet.dimord = '{chan}_spike_lfpchan_freq';
    %     morlet.trialtime = trialtimevec; 
    % 
    %     convol_nan = isnan(stsConvol.fourierspctrm{1,1});
    %     morlet_nan = morlet; 
    %     morlet_nan.fourierspctrm{1,1}(convol_nan) = NaN; 
    %     statSts2           = ft_spiketriggeredspectrum_stat(cfg,morlet);
    %     statSts3           = ft_spiketriggeredspectrum_stat(cfg,morlet_nan);

    %     subplot(1, 3, freq_lock);
    %     % plot the results
    %     hold on
    %     plot(statSts.freq,statSts.ppc1','r');
    %     %plot(statSts2.freq,statSts2.ppc0','g');
    %     plot(statSts3.freq,statSts3.ppc1','b');
    %     title(sprintf('%d Hz', lfp_freq(freq_lock)));
    %     xlabel('frequency, Hz');
    %     ylabel('ppc value');
    %     legend({'FT' 'Morlet'});
    %     set(gca,'XScale','log', 'XTick', [10 25 60]);
    %     xlim([min(keys.LFP.frequencies), max(keys.LFP.frequencies)]);
    
        %time_cmp(numtr, :) = [num_trials(numtr), ft_time, morlet_time, ft_time/morlet_time];
        
        time_cmp{1}(numtr, iter) = ft_time; % ft
        time_cmp{2}(numtr, iter) = morlet_time; % morlet
        time_cmp{3}(numtr, iter) = ft_time/morlet_time; %x times faster
        iter = iter+1;
    end
    
end

% calculate averages here
time_cmp{2, 1}(:, 1) = mean(time_cmp{1,1},2);
time_cmp{2, 1}(:, 2) = std(time_cmp{1,1},0,2);

time_cmp{2, 2}(:, 1) = mean(time_cmp{1,2},2);
time_cmp{2, 2}(:, 2) = std(time_cmp{1,2},0,2);

time_cmp{2, 3}(:, 1) = mean(time_cmp{1,3},2);
time_cmp{2, 3}(:, 2) = std(time_cmp{1,3},0,2);

pval = zeros(1,6);
for trnum = 1:6,
    [h,p,ci,stats] = ttest(time_cmp{1,1}(trnum, :), time_cmp{1,2}(trnum, :));
    pval(trnum) = p;
end
pval

save('ft_vs_morlet.mat', 'time_cmp');
% set(figure(h1), 'Paperunits','centimeters','PaperSize', wanted_papersize,'PaperPositionMode', 'manual','PaperPosition', [0 0 wanted_papersize])
% mtit(figure(h1),  [FigName_Long1 ], 'xoff', 0, 'yoff', 0.05, 'color', [0 0 0], 'fontsize', 14,'Interpreter', 'none');
% export_fig(figure(h1), [keys.path_to_save filesep FigName_short1], '-pdf','-transparent') % pdf by run