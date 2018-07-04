function bbedges = bayesianBlockBinEdges(data)
% Bayesian blocks
% This is a naive translation into Matlab
% from the Julia Implementation by Michael P.H. Stumpf and T. Chan
% Based on the Python Code of Jake Vanderplas.

% This version implements Bayesian blocks for histograms, where the
% data are sorted, then treated as event data (see Scargle 2012).
% Defaults as suggested in Scargle 2012 are used.

% References:
% Scargle 2012: http://adsabs.harvard.edu/abs/2012arXiv1207.5578S
% Python implementation: https://github.com/astroML/astroML/blob/master/astroML/density_estimation/bayesian_blocks.py
% Julia implementation: https://github.com/sisl/Discretizers.jl/blob/master/src/disc_bayesianblocks.jl
%
% INPUTS:
% data: should be a vector of values that you want to histogram

unique_data = unique(data); % sorted unique data
n_unique = length(unique_data);

bbedges = zeros(n_unique+1,1);
bbedges(1) = unique_data(1);
bbedges(2:end-1) = 0.5*(unique_data(1:end-1) + unique_data(2:end));
bbedges(end) = unique_data(end);
block_length = unique_data(end) - bbedges;

nn_vec = histcounts(data,bbedges);

count_vec = zeros(n_unique,1);
bestFit = zeros(n_unique,1);
lastInd = zeros(n_unique,1);

for kk = 1:n_unique
    widths = block_length(1:kk) - block_length(kk+1);
    count_vec(1:kk) = count_vec(1:kk) + nn_vec(kk);
    
    % Fitness function (eq. 19 from Scargle 2012)
    fit_vec = count_vec(1:kk).*log(count_vec(1:kk)./widths);
    % Prior (eq. 21 from Scargle 2012)
    fit_vec = fit_vec - (4 - log(73.53*0.05*kk^-0.478));
    fit_vec(2:end) = fit_vec(2:end) + bestFit(1:kk-1);
    
    [bestFit(kk), lastInd(kk)] = max(fit_vec);
end

change_points = zeros(n_unique,1);
ind_cp = n_unique+1;
indCurrent = n_unique+1;
while indCurrent >= 1
    ind_cp = ind_cp - 1;
    change_points(ind_cp) = indCurrent;
    if indCurrent==1
        break
    end
    indCurrent = lastInd(indCurrent - 1);
end
change_points = change_points(ind_cp:end);
bbedges = bbedges(change_points);

end