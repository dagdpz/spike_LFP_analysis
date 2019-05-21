function [ scores ] = MS_perm_test(data, keys, PPC_method)
%MS_T_TESTING Cluster-based corrected permutation t-test 
%   takes 2 vectors of values as a first variable
%   keys - for getting the vector of frequencies (time-points)
%   PPC_methods - lol, it seems like i am not even using it here

num_perm = 1000;
cr_value = 1.96; %ttest critical value for 0.05 p
sing_lev = 0.05;
clear scores

% frequencies as columns, trials as rows
pre_data = vertcat(data{1}(:, 1).ppc1);
post_data = vertcat(data{2}(:, 1).ppc1);

% Calculating p-values of the real data for defining clusters
for fr=1:numel(keys.LFP.frequencies),
    [h,p,ci,stats] = ttest(pre_data(:, fr), post_data(:, fr)); % h=1 - reject the null hypothesis
    scores.data(fr).p_val = p;
    scores.data(fr).tstat = stats.tstat;
end    
%plot(keys.LFP.frequencies, [scores.data.tstat]); 
%hline([1.96, -1.96]) 

% Identifying both positive and negative clusters
% Single values are kept as a cluster as well! (you can always exclude them
% later malually)
pos_tmp_ind = find([scores.data.tstat] > cr_value);
pos_difvec = diff(pos_tmp_ind);
pos_clusters = MS_find_clusters(pos_tmp_ind, pos_difvec);
    
neg_tmp_ind = find([scores.data.tstat] < -cr_value);
neg_difvec = diff(neg_tmp_ind);
neg_clusters = MS_find_clusters(neg_tmp_ind, neg_difvec); 

scores.clusters = [pos_clusters neg_clusters];

% Quit if there are no clusters in the real data
if isempty(scores.clusters),
    scores = [];
    return
end

% Calculating summed absolute t-test value for the real data within clusters
for clus = 1:numel(scores.clusters),
    if numel(scores.clusters{1, clus}) == 1,
        scores.clus_ttest(1, clus) = {abs(scores.data(1, scores.clusters{1, clus}).tstat)};
    else
        scores.clus_ttest(1, clus) = {sum(abs([scores.data(1, scores.clusters{1, clus}(1, 1):scores.clusters{1, clus}(1, 2)).tstat]))};
    end
end

% Finally I can permute the data!
scores.perm_ttest = zeros(num_perm, numel(scores.clusters));
data_full = vertcat(pre_data, post_data);
for prm=1:num_perm,
    perm_ind = randperm(numel(data_full(:,1)));
    
    for clus = 1:numel(scores.clusters),
        if numel(scores.clusters{1, clus}) == 1,
            [h,p,ci,stats] = ttest(data_full(perm_ind(1:numel(data_full(:, 1))/2), scores.clusters{clus}), ...
                data_full(perm_ind(numel(data_full(:, 1))/2+1:numel(data_full(:, 1))), scores.clusters{clus}));
            scores.perm_ttest(prm, clus) = abs(stats.tstat);
        else
            tsum = 0;
            for frv = scores.clusters{clus}(1,1):scores.clusters{clus}(1,2),
                [h,p,ci,stats] = ttest(data_full(perm_ind(1:numel(data_full(:, 1))/2), frv), ...
                    data_full(perm_ind(numel(data_full(:, 1))/2+1:numel(data_full(:, 1))), frv));
                tsum = tsum + abs(stats.tstat);
            end
            scores.perm_ttest(prm, clus) = abs(tsum);
        end
    end
end
scores.perm_dist = max(scores.perm_ttest, [], 2);
%hist(max(scores.perm_ttest, [], 2));

% Assign significance
for cl = 1:numel(scores.clus_ttest)
    tmpnum = sum(scores.clus_ttest{cl} > scores.perm_dist);
    pval = 1-tmpnum/num_perm;
    scores.out.pval(cl) = {pval};
    scores.out.h(cl) = {pval < sing_lev};
end


end