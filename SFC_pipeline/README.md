# Spike-Field Coherence analysis of the macaque parietal cortex activity within and across hemispheres  

## Main scripts  
**ph_run_LFP_analysis.m** - the function for ppc calculation (part of 'Lukas' DAG ephys pipeline)

Data in Y:\Projects\PPC_pulv_eye_hand\ephys\MIP_dPul_inj_working_post_sfn_SF

For each session, separate mat files for spikes and LFPs, 
output is spike_field_hands_inactivation_ppc1.mat

**MS_plot_PPC.m** - plot ppc before and after inactivation for each epoch and condition  
**MS_plot_ppc_perunit.m** - plot ppc before and after inactivation for each epoch and condition PER UNIT (averaged across LFP channels)  
**MS_scatter_per_unit.m** - plot *scatter plots* (one dot - one unit) before and after inactivation to see the overall change of tuning  

## Additional functions (called from MS_plot_PPC)  
**MS_perm_test.m** - Cluster-adjusted permutation t-test (as described in Maris & Oostenveld, 2007)  
**MS_find_clusters.m** - accessory function for finding clusters in permutation test  

plus

**MS_FDR_test.m** - False Discovery Rate (FDR) correction for multiple comparisons (NOT USED IN MS thesis) 


### Example input  

**ph_run_LFP_analysis.m** C:\Users\mslashcheva\Dropbox\DAG\DAG_toolbox\spike_analysis as following:  

**ph_run_LFP_analysis('PPC_pulv_eye_hand',{'MIP_dPul_inj_working_post_sfn'})**

Run this in a script or in a command prompt:

Loc = 'Y:\Projects\PPC_pulv_eye_hand\ephys\MIP_dPul_inj_working_post_sfn_SF\spike_field_hands_inactivation_ppc1.mat'; %MIP  
PPC_method='ppc1';  
hemis = {'MIP_L', 'MIP_R', 'MIP_L_MIP_R', 'MIP_R_MIP_L'};  

for hemi = hemis,  
    disp(hemi)  
    hemi = hemi{:};  
    MS_plot_PPC_v1(Loc, PPC_method, hemi)  
end  
