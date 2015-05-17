%% linear_mixing.m
% creates a linear combination of a number of Gaussians with random parameters drawn from a larger library, and
% then inverts to recover the coefficients of each contribution.
close all
clear all

addpath('C:\Documents\Work\Code\MATLAB\utilities')


%% define parameters
    % number of species in library
    library_n = 50;
    % number of species in the sample
    sample_n = 10;

% create x-axis
    kmin = 1000;
    kmax = 4000;
    k_n = 200;
    K0 = linspace(kmin, kmax, k_n);
    
% Define the spectral response of each species, and assemble into a matrix.
     cmax = 0.95*kmax;
     cmin = 1.05*kmin;
     wmax = 0.1*(kmax - kmin);
     wmin = 0.01*(kmax - kmin);
     hmax = 1;
     hmin = 0.1;

     centers = (cmax - cmin)*rand([1 library_n]) + cmin;
     widths = (wmax - wmin)*rand([1 library_n]) + wmin;
     heights = (hmax - hmin)*rand([1 library_n]) + hmin;
     
    
     
for m = 1:library_n;
    spectra_lib(m,:) = heights(m).*exp(-(K0 - centers(m)).^2/(2*widths(m).^2));
end

% plot all spectra in the library
plot_single(K0,spectra_lib',kmin,kmax,'Wavenumber [cm^{-1}]','Transmittance')
%plot(1-spectra_lib')

%% choose nspec of these as the species in the sample
p = randperm(library_n);
index_in = sort(p(1:sample_n));
spectra_in = spectra_lib(index_in,:)';

% plot all spectra in the sample
plot_single(K0,spectra_in,kmin,kmax,'Wavenumber [cm^{-1}]','Transmittance')
%plot(1 - spectra_in)


% plot composite spectrum
min_dens = 0.01;
densities = rand(1,sample_n)';
densities = max(densities/(sum(densities)),min_dens);
composite_spectrum = spectra_in*densities;

dens_matrix = repmat(densities,1,k_n)';
plot_spectsum(K0,cat(2,dens_matrix.*spectra_in,composite_spectrum),kmin,kmax,'Wavenumber [cm^{-1}]','Transmittance')
%plot_single(K0,1 - dens_matrix.*spectra_in,kmin,kmax,'Wavenumber [cm^{-1}]','Transmittance')


%figure,
%plot(1-composite_spectrum)

% add noise
snr = 1000;
noise = normrnd(0,1/snr,[1,k_n])';
% signal = composite_spectrum + noise;
signal = composite_spectrum;

% plot(1-abs(signal))

% invert
densities_out = pinv(spectra_lib')*signal;

% identify species present
thresh = 0.005;
found_indices = densities_out > thresh;
index1 = 1:library_n;
found_species = index1(found_indices);

% display results
% index_in
% found_species

% plot 
plot_single(K0,composite_spectrum,kmin,kmax,'Wavenumber [cm^{-1}]','Transmittance')
plot_library(K0,repmat(1:library_n,k_n,1),spectra_lib)



% [FileName,PathName] = uiputfile('*.*','Save Absorbance Spectrum and .eps Figure')
% print2eps(strcat(PathName,FileName,'.eps'),handle1)