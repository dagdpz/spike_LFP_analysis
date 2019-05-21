# Spike-field coherence simulation
The collection of scripts for the simulation of spike-field coherence

## Main scripts

**MS_SFCsim_coh_comparison.m** - main script for SFC simulation and comparison of Spectral coherence, calculated with Chronux and FieldTrip, and Pairwise Phase consistency (ppc)  
**MS_SFCsim_coh_comparison_bursting.m** - same as previous, but with locked spikes firing in bursts  

**MS_SFCsim_FT_vs_Morlets.m** - comparison of Fourier transform and Complex Morlet Wavelet calculation of phase around each spike and the resulting ppc  
**MS_SFCsim_num_of_observations.m** - Simulation of SFC with different number of pairs of spike-field phases  


## Additional scripts

**MS_SFCsim_basic.m** - the very basic implementation of simulation  
**MS_SFCsim_bursting.m** - the implementation of bursting spikes  

**MS_test_FT_simulate_spike_train.m** - accessory function for the generation of Poissanian spike train  
**test_FT_fieldtrip2chronux.m** - accessory function for converting data from FieldTrip format to Chronux  

**MS_plotTrial.m** - function for vilualizing the simulation  
**MS_plotRaster.m** - accessory function for raster plot  

**MS_test_FT_Morlet_convol_basic.m** - the very basic implementation of Complex Morlet Wavelet convolution  

