function [fr,edges,spike_train] = calc_ratemod(cl_evix,ev_ix,sr,binsize, varargin);
% Calc rate modulation


gauswin = 25;
% binsize = 0.01;

    
evL = length(ev_ix)/sr;
spks = (cl_evix-ev_ix(1))./sr;
edges = [0:binsize:evL];


[fr,spike_train]  = deal(zeros(1,length(edges)-1));
if isempty(cl_evix); return; end;


spike_train = histcounts(spks,edges);
gauss_kernel = gausswin(gauswin);                   % n-sample Gaussian kernel
gauss_kernel = gauss_kernel/sum(gauss_kernel);                    % normalize

% % Convert sigma to samples
% sigma_samples = sigma_sec / binsize;

% % Set kernel to cover ±3 sigma
% kernel_half_width = round(3 * sigma_samples);
% kernel_len = 2 * kernel_half_width + 1;
% gauss_kernel = gausswin(kernel_len);                   % n-sample Gaussian kernel



fr = conv(spike_train, gauss_kernel, 'same')/ binsize ; % smoothed firing rate (Hz)


end