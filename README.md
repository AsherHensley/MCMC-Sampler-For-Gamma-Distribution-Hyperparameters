# MCMCSamplerForGammaHyperparameters
Sample posterior of gamma distribution Hyperparameters with Metropolis Hastings algorithm.

Call main.m to run on simulated data

%METROPOLISHASTINGS(X,N) Sample gamma hyperparameters w/ Metropolis-Hastings
%   [A,B] = SAMPLEGAMMAPARM(X,N) returns N samples from the posterior of
%   the gamma distribution's shape and rate parameters A and B given a set 
%   of observations X using the Metropolis-Hastings algorithm. The proposal
%   distribution used is an isotropic truncated Gaussian distribution
%   centered on the previous [A,B] sample. Truncation is used to make sure 
%   no negative values of [A,B] are proposed. 
%
%   The observations are assumed to be drawn from the following model:
%
%   [A,B] ~ p(A,B)
%       X ~ Gamma(X|A,B)
%
%   The program then returns N samples from the posterior:
%
%   [A,B] ~ p(A,B|X)
%
%   [A,B,REJECT] = SAMPLEGAMMAPARM(X,N) returns the REJECT mask and
%   indicating which samples of A and B were rejected by the 
%   Metropolis-Hastings update.
%
%   [A,B,REJECT] = SAMPLEGAMMAPARM(...,'PARAM',VALUE) allows the use
%   of optional input parameters:
%       'sigma' - Use custom sigma for the proposal distribution [default = 0.1]. 
%       'icond' - Set the initial state of the Markov chain (default = [0;0]).
