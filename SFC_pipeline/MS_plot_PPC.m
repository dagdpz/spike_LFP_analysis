function MS_plot_PPC( spike_field_location, PPC_method, hemi )
%	Inputs: SFfile_location (e.g. spike_field_hands_inactivation_ppc1.mat), PPC_method, hemi

%   PPC_method is 'ppc0' or 'ppc1'
%   hemi is MIP_L, MIP_R, MIP_R_MIP_L, MIP_L_MIP_R

% Output: plots of ppc values before and after inactivaion for each
% condition and epoch separately

load(spike_field_location);

plot_from_perturbation_effect = 1;
perturbation_group = 2;
calculate_ttest = 1;
spike_pairs = 2500; 
bands_v = [4, 8, 12, 30, 100];
SF_combinations = [spike_field.unit];

% Removing units not included in SFN analysis 
% SFN list has 78 cells, here we have 84, 
remove_units = 0; 
if remove_units == 1,
    pairs_sep_cells = [{SF_combinations.unit_ID}]';
    cell_list_sfn = importdata('Y:\Projects\PPC_pulv_eye_hand\ephys\MIP_dPul_inj_working_post_sfn\cell_list.txt');
    disp(sprintf('the number of units in the current data is %d, in the SFN data - %d', numel(unique(pairs_sep_cells)), numel(cell_list_sfn))); 
    disp('Units to be excluded');
    unique(pairs_sep_cells(ismember(pairs_sep_cells, cell_list_sfn) == 0))
    SF_combinations(find((ismember(pairs_sep_cells, cell_list_sfn) == 0))) = [];
end

siteIDs_tmp = arrayfun(@(x) repmat({x.site_ID}, numel(x.unit),1), spike_field, 'uniformoutput',false);
LFP_siteIDs = vertcat(siteIDs_tmp{:});

hemi_tmp = arrayfun(@(x) repmat({x.target}, numel(x.unit),1), spike_field, 'uniformoutput',false);
LFP_hemiIDs = (vertcat(hemi_tmp{:}));

if remove_units == 1;
    LFP_siteIDs(find((ismember(pairs_sep_cells, cell_list_sfn) == 0))) = [];
    LFP_hemiIDs(find((ismember(pairs_sep_cells, cell_list_sfn) == 0))) = [];
end


% getting info about the same site unit and LFP, 1-same, 0-different
same_site_ind = strcmp(LFP_siteIDs, {SF_combinations.site_ID}');

% getting info about hemisphere - ACROSS - LFP and spikes are in different
% hemispheres
for ind=1:numel(LFP_hemiIDs),
    if strcmp(LFP_hemiIDs(ind), SF_combinations(ind).target),
        same_hemi_ind_t(ind) = LFP_hemiIDs(ind);
    else
        same_hemi_ind_t(ind)= {[SF_combinations(ind).target, '_', LFP_hemiIDs{ind}]}; %{'ACROSS'};
    end
end
same_hemi_ind = same_hemi_ind_t';

% choose hemisphere AND remove same-site pairs
subset_ind = strcmp(same_hemi_ind, {hemi}) & ~same_site_ind;

% subset of all the pairs from indicated hemi and different sites
SF_difchan_hemi = SF_combinations(subset_ind);
%(BG) Get current hemisphere and load avg unit struct from
%perturbation table
current_target = SF_difchan_hemi.target;
load (['Y:\Projects\PPC_pulv_eye_hand\ephys\MIP_dPul_inj_w_o_20171012\perturbation_table\Linus_' current_target '_Ddre_han_avg_units_struct']);
LFP_siteID_difchan_hemi = LFP_siteIDs(subset_ind); %LFP channels

% checking the number of unique units (for debug)
% [b,m,n] = unique([{SF_difchan_hemi.unit_ID}]');
% single_rate = [SF_difchan_hemi(m).Single_rating];
% sum(single_rate==1) %SU
% numel(single_rate) - sum(single_rate==1) %MU

%% Conditions to compare keys
keys.conditions_to_compare{1}(1).reach_hand=1;
keys.conditions_to_compare{1}(1).hemifield=-1;
keys.conditions_to_compare{1}(1).choice=0;
keys.conditions_to_compare{1}(1).perturbation=0;
keys.conditions_to_compare{1}(1).color='b';
keys.conditions_to_compare{1}(1).title='LHLS';

keys.conditions_to_compare{1}(2).reach_hand=1;
keys.conditions_to_compare{1}(2).hemifield=-1;
keys.conditions_to_compare{1}(2).choice=0;
keys.conditions_to_compare{1}(2).perturbation=1;
keys.conditions_to_compare{1}(2).color='r';
keys.conditions_to_compare{1}(2).title='LHLS';


keys.conditions_to_compare{2}(1).reach_hand=1;
keys.conditions_to_compare{2}(1).hemifield=1;
keys.conditions_to_compare{2}(1).choice=0;
keys.conditions_to_compare{2}(1).perturbation=0;
keys.conditions_to_compare{2}(1).color='b';
keys.conditions_to_compare{2}(1).title='LHRS';

keys.conditions_to_compare{2}(2).reach_hand=1;
keys.conditions_to_compare{2}(2).hemifield=1;
keys.conditions_to_compare{2}(2).choice=0;
keys.conditions_to_compare{2}(2).perturbation=1;
keys.conditions_to_compare{2}(2).color='r';
keys.conditions_to_compare{2}(2).title='LHRS';


keys.conditions_to_compare{3}(1).reach_hand=2;
keys.conditions_to_compare{3}(1).hemifield=-1;
keys.conditions_to_compare{3}(1).choice=0;
keys.conditions_to_compare{3}(1).perturbation=0;
keys.conditions_to_compare{3}(1).color='b';
keys.conditions_to_compare{3}(1).title='RHLS';

keys.conditions_to_compare{3}(2).reach_hand=2;
keys.conditions_to_compare{3}(2).hemifield=-1;
keys.conditions_to_compare{3}(2).choice=0;
keys.conditions_to_compare{3}(2).perturbation=1;
keys.conditions_to_compare{3}(2).color='r';
keys.conditions_to_compare{3}(2).title='RHLS';


keys.conditions_to_compare{4}(1).reach_hand=2;
keys.conditions_to_compare{4}(1).hemifield=1;
keys.conditions_to_compare{4}(1).choice=0;
keys.conditions_to_compare{4}(1).perturbation=0;
keys.conditions_to_compare{4}(1).color='b';
keys.conditions_to_compare{4}(1).title='RHRS';

keys.conditions_to_compare{4}(2).reach_hand=2;
keys.conditions_to_compare{4}(2).hemifield=1;
keys.conditions_to_compare{4}(2).choice=0;
keys.conditions_to_compare{4}(2).perturbation=1;
keys.conditions_to_compare{4}(2).color='r';
keys.conditions_to_compare{4}(2).title='RHRS';

epochs = [keys.EPOCHS_PER_TYPE{4}(:,1)'];


clear all_cond
Parameters={'type','effector','reach_hand','choice','perturbation','hemifield'};
for par=1:numel(Parameters),
    SFidx.(Parameters{par})=arrayfun(@(x) isfield(x.per_condition,Parameters{par}) && numel(unique([x.per_condition.(Parameters{par})]))>1,SF_difchan_hemi);
end

limit_pairs_by_condition={'perturbation'};
valid_SF_combinations=true(size(SF_difchan_hemi));
for L=1:numel(limit_pairs_by_condition)
    valid_SF_combinations = valid_SF_combinations &   SFidx.(limit_pairs_by_condition{L});
end

% Plotting population ppc
n_unique_units=numel(unique({SF_difchan_hemi(valid_SF_combinations).unit_ID}));
n_unique_sites=numel(unique({SF_difchan_hemi(valid_SF_combinations).site_ID}));
disp(['unique_sites:' num2str(n_unique_sites) ' ,unique_units:' num2str(n_unique_units)]);


%% Finding the same set of SF pairs within each epoch

SF_difchan_hemi_tmp = SF_difchan_hemi(valid_SF_combinations);
SF_same_set = zeros(sum(valid_SF_combinations), numel(epochs));
for un=1:numel(SF_difchan_hemi_tmp),
    un_tmp = [SF_difchan_hemi_tmp(un).per_condition];
    un_tmp = un_tmp(logical(~[un_tmp.choice]));
    num_p = zeros(numel(un_tmp), numel(epochs));
 
    for cnd=1:numel(un_tmp),
        ep_tmp = [un_tmp(cnd).per_epoch];
        num_p(cnd, :) = [ep_tmp.n_pairs];
    end
    SF_same_set(un, :) = all(num_p > spike_pairs, 1);
end
% sum(SF_same_set)

%
keys.path_to_save = sprintf('%s\\%s\\SFC_population_%s', [keys.drive, keys.basepath_to_save], keys.project_version, PPC_method);

if exist(keys.path_to_save) == 0;
    mkdir(keys.path_to_save);
end

if plot_from_perturbation_effect
 FigName_short2= sprintf('%s_%s_%d_exc_same_pairs_STAT_%d', hemi,PPC_method, spike_pairs,perturbation_group); %name of the file
FigName_Long2=sprintf('%s_%s, %d/%d SF pairs, %d spike pairs, perturbation effect %d', ...
    hemi, PPC_method, sum(valid_SF_combinations), numel(SF_difchan_hemi), spike_pairs, perturbation_group ); % overall title of the figure   
else
FigName_short2= sprintf('%s_%s_%d_exc_same_pairs_STAT', hemi,PPC_method, spike_pairs); %name of the file
FigName_Long2=sprintf('%s_%s, %d/%d SF pairs, %d spike pairs', ...
    hemi, PPC_method, sum(valid_SF_combinations), numel(SF_difchan_hemi), spike_pairs); % overall title of the figure
end
set(0,'DefaultTextInterpreter','none');
wanted_papersize= [40 25];
h1=figure('outerposition',[0 0 1900 1200],'name', FigName_Long2);
figure(h1)

clear H
clear perm_scores
perm_scores.hemi = hemi;
if ~calculate_ttest,
    load(sprintf('%s\\perm_scores_%s.mat', keys.path_to_save, hemi))
end
%% PLOT

for ep=1:numel(epochs)
    y_max = 0;
    y_min = 0;
    y_bar = 0;
    max_cond = 1;
    unique_un = numel(unique({SF_difchan_hemi_tmp(logical(SF_same_set(:, ep))).unit_ID}));
    
    
    if plot_from_perturbation_effect 
    for un=1:numel(SF_difchan_hemi_tmp)
       
            if any(strfind(current_target,'_L'))
                un_ind_tmp = find(strcmp(avg_unit_epoch_MIP_L.unit_ID, SF_difchan_hemi_tmp(1,un).unit_ID));    
                pert_group_tmp(un,1) = avg_unit_epoch_MIP_L.(epochs{ep})(un_ind_tmp,1) == perturbation_group;
            else
                un_ind_tmp = find(strcmp(avg_unit_epoch_MIP_R.unit_ID, SF_difchan_hemi_tmp(1,un).unit_ID));    
                pert_group_tmp(un,1) = avg_unit_epoch_MIP_R.(epochs{ep})(un_ind_tmp,1) == perturbation_group; 
            end   
    end
    SF_comb_per_cond = [SF_difchan_hemi_tmp((logical(SF_same_set(:, ep)) & logical(pert_group_tmp))).per_condition];
    unique_un = numel(unique({SF_difchan_hemi_tmp((logical(SF_same_set(:, ep)) & logical(pert_group_tmp))).unit_ID}));
    else
    
    SF_comb_per_cond = [SF_difchan_hemi_tmp(logical(SF_same_set(:, ep))).per_condition];
    end
    
%     %Take only unique units present in the pairs and avg struct and their
%     %avg_unit response from pert_table (BG)
%     unit_ID_unique_pairs = unique({SF_difchan_hemi_tmp.unit_ID});
%     common_units_1 = struct('unit_ID', [], 'avg_unit_response',[]);
%     common_units_2 = struct('unit_ID', [], 'avg_unit_response',[]);
%     common_units_0 =struct('unit_ID', [], 'avg_unit_response',[]);
%     counter0=1; 
%     counter1=1;
%     counter2=1;
%     
%     for unit_IDi=1:numel (unit_ID_unique_pairs)
%         if any(strfind(current_target,'_L'))
%             unit_ID_unique_avg = unique ({avg_unit_epoch_MIP_L.unit_IDs_L});
%             if cellfun(@strcmp, unit_ID_unique_pairs(unit_IDi), unit_ID_unique_avg (unit_IDi)) 
%                 if avg_unit_epoch_MIP_L(unit_IDi).(epochs{ep})==0
%                     common_units_0(counter0).unit_ID  = unit_ID_unique_pairs(unit_IDi);
%                     common_units_0(counter0).avg_unit_response = 0;
%                     counter0=counter0+1;
%                 elseif avg_unit_epoch_MIP_L(unit_IDi).(epochs{ep})==1
%                     common_units_1(counter1).unit_ID  = unit_ID_unique_pairs(unit_IDi);
%                     common_units_1(counter1).avg_unit_response = 1;
%                     counter1=counter1+1;
%                 elseif avg_unit_epoch_MIP_L(unit_IDi).(epochs{ep})==2
%                     common_units_2(counter2).unit_ID  = unit_ID_unique_pairs(unit_IDi);
%                     common_units_2(counter2).avg_unit_response = 2;
%                     counter2=counter2+1;
%                 end
%             end
%         else
%             unit_ID_unique_avg = unique ({avg_unit_epoch_MIP_R.unit_IDs_R});
%             if cellfun(@strcmp, unit_ID_unique_pairs(unit_IDi), unit_ID_unique_avg (unit_IDi)) 
%                if avg_unit_epoch_MIP_R(unit_IDi).(epochs{ep})==0
%                     common_units_0(unit_IDi).unit_ID  = unit_ID_unique_pairs(unit_IDi);
%                     common_units_0(unit_IDi).avg_unit_response = 0;
%                 elseif avg_unit_epoch_MIP_R(unit_IDi).(epochs{ep})==1
%                     common_units_1(unit_IDi).unit_ID  = unit_ID_unique_pairs(unit_IDi);
%                     common_units_1(unit_IDi).avg_unit_response = 1; 
%                 elseif avg_unit_epoch_MIP_R(unit_IDi).(epochs{ep})==2
%                     common_units_2(unit_IDi).unit_ID  = unit_ID_unique_pairs(unit_IDi);
%                     common_units_2(unit_IDi).avg_unit_response = 2;
%                 end
%             end
%         end
%     end

    clear all_cond
    if isempty(SF_comb_per_cond)
       continue 
    end
    for par=1:numel(Parameters),
        all_cond(:, par) = [SF_comb_per_cond.(Parameters{par})];
        paridx.(Parameters{par}) = par;
    end

    for row=1:numel(keys.conditions_to_compare)
        condition_fieldnames=fieldnames(keys.conditions_to_compare{row});
        condition_fieldnames=condition_fieldnames(ismember(condition_fieldnames,Parameters));
        
        % here i collect pre-and post of one condition
        if calculate_ttest,
        clear perm_test
            for con=1:numel(keys.conditions_to_compare{row}) % this is ugly, but the same loop again
                current_index=true(size(all_cond,1),1);
                for FN=condition_fieldnames'
                    current_index=current_index & all_cond(:, paridx.(FN{:}))==keys.conditions_to_compare{row}(con).(FN{:});
                end
                perm_test_tmp = vertcat(SF_comb_per_cond(current_index).per_epoch);
                perm_test(con) = {perm_test_tmp(:, ep)}; % not very smart
            end
            % Here will be the function with the permutation testing
            clear t_scores
            t_scores = MS_perm_test(perm_test, keys, PPC_method);
            perm_scores.score(row, ep) = {t_scores};
        else
            t_scores = perm_scores.score{row, ep};
        end
        
        for con=1:numel(keys.conditions_to_compare{row}) % this is ugly, but the same loop again
            current_index=true(size(all_cond,1),1);
            for FN=condition_fieldnames'
                current_index=current_index & all_cond(:, paridx.(FN{:}))==keys.conditions_to_compare{row}(con).(FN{:});
            end
            SF_comb_per_cond_per_epoch = vertcat(SF_comb_per_cond(current_index).per_epoch);
            SF_comb_per_cond_per_epoch = SF_comb_per_cond_per_epoch(:, ep); % not very smart

            ax(row) = subplot(numel(keys.conditions_to_compare), numel(epochs), numel(epochs)*(row-1)+ep); hold on;
            
            %temptemp=arrayfun(@(x) numel(x.(PPC_method)),SF_comb_per_cond_per_epoch(:,ep))>1;
            PPC=vertcat(SF_comb_per_cond_per_epoch.(PPC_method));
            PPC(isinf(PPC))=NaN;
            
            if max(mean(PPC, 1)) > y_max, 
                y_max = max(mean(PPC, 1) + 3*max(sem(PPC,1))); 
                y_bar = max(mean(PPC, 1) + 1*max(sem(PPC,1))); 
                max_cond = row; 
            end
            if min(mean(PPC, 1)) < y_min, y_min = min(mean(PPC, 1) - max(sem(PPC,1))); end

            H(ep) = shadedErrorBar(keys.LFP.frequencies, mean(PPC, 1), sem(PPC,1), keys.conditions_to_compare{row}(con).color,1);
            y_lim=get(gca,'ylim');
            set(gca,'XScale','log') 
            
            % text(min(keys.LFP.frequencies),y_lim(2)-con*diff(y_lim)/10,['Excluded:' num2str(sum(all(isnan(PPC),2)) + sum(~temptemp))],'color', keys.conditions_to_compare{row}(con).color);
            if row==1 && con==2
                title(sprintf('%s, %d ms \nSF pairs: %d; Units: %d \n', epochs{ep}, round((keys.EPOCHS_PER_TYPE{4}{ep, 4} - keys.EPOCHS_PER_TYPE{4}{ep, 3})*1000),...
                    numel(PPC(:, 1)), unique_un));
%                 text(min(keys.LFP.frequencies),y_lim(2)-con*diff(y_lim)/10,...
%                     ['SF_pairs:' num2str(numel(PPC(:, 1))) newline 'Units:' num2str(unique_un)], 'Color', 'r');
            end
            
            if con==1 && ep==1
                ylabel([keys.conditions_to_compare{row}(con).title]);
            end
            
            if con == 2,

                set(gca,'XTick',[bands_v]); %round(10.^(0.6:0.2:2))
%                 set(gca,'YGrid', 'off');
                if ~isempty(t_scores),
                    %plot_clust = t_scores.clusters;
                    plot_clust = t_scores.clusters([t_scores.out.h{1,:}]);
                    for cl=1:numel(plot_clust), %select clusters here - ([t_scores.out.h{1,:}])
                        if numel(plot_clust{cl}) == 1;
                            %text(keys.LFP.frequencies(1, plot_clust{cl}), y_bar, ['*'], 'Color', [0.7 0.7 0.7], 'FontSize',14);
                            continue
                        else
                            line([keys.LFP.frequencies(1, plot_clust{cl}(1,1)) keys.LFP.frequencies(1, plot_clust{cl}(1,2))], ...
                                [y_bar y_bar], 'Color', [0.7 0.7 0.7], 'LineWidth', 4); 
                        end
                    end
                end
            end
            xlim([min(keys.LFP.frequencies), max(keys.LFP.frequencies)]);
        end
    end
    set(ax(1:4), 'ylim', [y_min y_max]);
    for sp = 1:row,
        subplot(ax(sp));
        z = vline([bands_v]); %, 'r:'
        set(z, 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5, 'LineStyle', ':');
    end
end

set(figure(h1), 'Paperunits','centimeters','PaperSize', wanted_papersize,'PaperPositionMode', 'manual','PaperPosition', [0 0 wanted_papersize])
mtit(figure(h1),  [FigName_Long2 ], 'xoff', 0, 'yoff', 0.05, 'color', [0 0 0], 'fontsize', 12,'Interpreter', 'none');
export_fig(figure(h1), [keys.path_to_save filesep FigName_short2], '-pdf','-transparent') % pdf by run

scores_name = [keys.path_to_save filesep 'perm_scores_' hemi '.mat'];
save(scores_name,'perm_scores');

close all

end

% Plotting part for the initial implementation
% for row=1:numel(keys.conditions_to_compare)
%     condition_fieldnames=fieldnames(keys.conditions_to_compare{row});
%     condition_fieldnames=condition_fieldnames(ismember(condition_fieldnames,Parameters));
%     
%     no_ppc_ind = zeros(sum(valid_SF_combinations), numel(epochs));
%     
%     for con=1:numel(keys.conditions_to_compare{row})
%         
%         current_index=true(size(all_cond,1),1);
%         for FN=condition_fieldnames'
%             current_index=current_index & all_cond(:, paridx.(FN{:}))==keys.conditions_to_compare{row}(con).(FN{:});
%         end
%         
%         SF_comb_per_cond_per_epoch = vertcat(SF_comb_per_cond(current_index).per_epoch);
%         check for NaN values and only one value of ppc - assigning 1;
%         SF_no_ppc = arrayfun(@(x) all(isnan(x.(PPC_method)))| numel(x.(PPC_method)) < 2 ...
%             | x.n_pairs<spike_pairs, SF_comb_per_cond_per_epoch);
%         no_ppc_ind = no_ppc_ind | SF_no_ppc; % indices of pairs for exclusion based on low FR
%     end
%     
%     This is for the statistical testing and permutations which runs as a
%     separate function
%     clear stat_test
%     for con=1:numel(keys.conditions_to_compare{row}) % this is ugly, but the same loop again
%         current_index=true(size(all_cond,1),1);
%         for FN=condition_fieldnames'
%             current_index=current_index & all_cond(:, paridx.(FN{:}))==keys.conditions_to_compare{row}(con).(FN{:});
%         end
%         SF_comb_per_cond_per_epoch = vertcat(SF_comb_per_cond(current_index).per_epoch);
%         stat_test(con) = {SF_comb_per_cond_per_epoch};
%     end    
%     t_scores = MS_t_testing(stat_test, no_ppc_ind, 1000, keys, epochs, PPC_method);
%     perm_scores.scores(row) = {t_scores};
%     
%     Plotting is actually happening here
%     for con=1:numel(keys.conditions_to_compare{row}) % this is ugly, but the same loop again
%         current_index=true(size(all_cond,1),1);
%         for FN=condition_fieldnames'
%             current_index=current_index & all_cond(:, paridx.(FN{:}))==keys.conditions_to_compare{row}(con).(FN{:});
%         end
%         SF_comb_per_cond_per_epoch = vertcat(SF_comb_per_cond(current_index).per_epoch);
%         
%         for ep=1:numel(epochs),
%             subplot(numel(keys.conditions_to_compare), numel(epochs), numel(epochs)*(row-1)+ep); hold on;
%             
%             temptemp=arrayfun(@(x) numel(x.(PPC_method)),SF_comb_per_cond_per_epoch(:,ep))>1;
%             PPC=vertcat(SF_comb_per_cond_per_epoch(~no_ppc_ind(:, ep),ep).(PPC_method));
%             PPC(isinf(PPC))=NaN;
%             nan_num = sum(all(isnan(PPC), 2));
%             
%             H(ep) = shadedErrorBar(keys.LFP.frequencies, mean(PPC, 1), sem(PPC,1), keys.conditions_to_compare{row}(con).color);
%             y_lim=get(gca,'ylim');
%             set(gca,'XScale','log') 
%             text(min(keys.LFP.frequencies),y_lim(2)-con*diff(y_lim)/10,['Excluded:' num2str(sum(all(isnan(PPC),2)) + sum(~temptemp))],'color', keys.conditions_to_compare{row}(con).color);
%             text(min(keys.LFP.frequencies),y_lim(2)-con*diff(y_lim)/10,['pairs:' num2str(numel(PPC(:, 1)))],'color', keys.conditions_to_compare{row}(con).color);
% 
%                         hold on;
%                         H(ep) = shadedErrorBar(5:100, mean(nopert_tmp, 1), sem(nopert_tmp), 'b');
%             if con==1
%                 xlabel('frequency');
%                 title(sprintf('%s, %d ms', epochs{ep}, round((keys.EPOCHS_PER_TYPE{4}{ep, 4} - keys.EPOCHS_PER_TYPE{4}{ep, 3})*1000) )); % ???????
%             end
%             
%             if con==2
%                 vline(bands_v); %vline([4, 8, 12, 30]);
%                 for band = 1:4,
%                     if t_scores(ep).out(band).signif == 1,
%                         text(bands(band), 0, sprintf('* %d',t_scores(ep).out(band).pvalue));
%                         text(bands_v(band), max(mean(PPC, 1)), sprintf('*'), 'FontSize',18);
%                     end
%                 end
%             end
%             
%             if con==1 && ep==1
%                 ylabel([keys.conditions_to_compare{row}(con).title ', N=' num2str(sum(current_index))]);
%             end
%             ylim([-0.001 0.015]);
%             xlim([min(keys.LFP.frequencies), max(keys.LFP.frequencies)]);
%         end
%     end
% end
% 
% set(figure(h1), 'Paperunits','centimeters','PaperSize', wanted_papersize,'PaperPositionMode', 'manual','PaperPosition', [0 0 wanted_papersize])
% mtit(figure(h1),  [FigName_Long2 ], 'xoff', 0, 'yoff', 0.05, 'color', [0 0 0], 'fontsize', 14,'Interpreter', 'none');
% export_fig(figure(h1), [keys.path_to_save filesep FigName_short2], '-pdf','-transparent') % pdf by run
% 
% scores_name = [keys.path_to_save filesep 'perm_scores_' hemi '.mat'];
% save(scores_name,'perm_scores');
% 
% close all
% 
% end


