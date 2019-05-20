function [ output_args ] = MS_plotTrial( which_trial, data, spikePhaseFreq, spikeRate, spikeLockProb, trllen, lockedSpikes, Spikes )
% Use this function to take a look at how simulated data looks like

ax1 = subplot(4, 1, 1);
plot(data.time{1}, data.trial{which_trial}(1, :)); % mixed signal
ylabel('LFP amplitude', 'fontsize',8);
title('Superimposed LFP signal');
set(gca,'xtick',[]);
box(ax1,'off')

ax2 = subplot(4,1,2);
plot(data.time{1}, data.trial{which_trial}(spikePhaseFreq+1, :)); % LFP channel with locked spikes
ylabel('LFP amplitude', 'fontsize',8);
title('LFP component with locked spikes (10Hz)');
set(gca,'xtick',[]);
box(ax2,'off')


ax3 = subplot(4,1,3);
MS_plotRaster(logical(lockedSpikes), data.time{1});
ylabel('Trial number', 'fontsize',8);
set(gca, 'ytick', [0 10 20 30]);
title('Locked spikes only (bursts)')
% title(sprintf('Locked; FR %d Hz, LockProb %.2f, mean locked fr %.2f, %d locked spikes', ... 
%     spikeRate,spikeLockProb,mean(sum(lockedSpikes,2)/trllen),sum(sum(lockedSpikes))));
set(gca,'xtick',[]);


ax4 = subplot(4,1,4);
MS_plotRaster(Spikes, data.time{1});
xlabel('Time, sec', 'fontsize',8);
ylabel('Trial number', 'fontsize',8);
set(gca, 'ytick', [0 10 20 30]);
title('All spikes (locked and random)')
% title(sprintf('All spikes; FR %d Hz, spikeLockProb %.2f, mean rate %.2f, %d spikes',...
%     spikeRate,spikeLockProb,mean(sum(Spikes,2)/trllen),sum(sum(Spikes))));

end

