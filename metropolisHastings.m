function [a,b,reject] = metropolisHastings(x,n,varargin)
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
%
%   Copyright 2016 Asher A Hensley
%   $Revision: 1.0 $  $Date: 2016/09/30 $
%
%   MIT License
%   
%   Permission is hereby granted, free of charge, to any person obtaining a 
%   copy of this software and associated documentation files (the 
%   "Software"), to deal in the Software without restriction, including 
%   without limitation the rights to use, copy, modify, merge, publish, 
%   distribute, sublicense, and/or sell copies of the Software, and to 
%   permit persons to whom the Software is furnished to do so, subject to 
%   the following conditions:
% 
%   The above copyright notice and this permission notice shall be included 
%   in all copies or substantial portions of the Software.
% 
%   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS 
%   OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
%   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
%   IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY 
%   CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
%   TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
%   SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

%Setup
X = zeros(2,n);
SIGMA = 0.1;
reject = zeros(1,n);
ICOND = [0;0];

%Optional Input Args
if nargin>2
    narg = length(varargin);
    for kk = 1:2:narg
        switch varargin{kk}
            case 'sigma'
                SIGMA = varargin{kk+1};   
            case 'icond'
                ICOND = varargin{kk+1}; 
            otherwise
                error('Unknown parameter')
        end
    end
end

%Set Initial State
X(:,1) = ICOND(:);

%Proposal Distribution
Q = @(z,mu,sig)mvnpdf(z,mu,sig^2*eye(2))/(mvncdf([0;0],inf(2,1),mu,sig^2*eye(2)));

%Target Distribution
p = prod(x);
q = sum(x);
r = length(x);
s = length(x);
P = @(z)p.^(z(1)-1).*exp(-z(2)*q)./(gamma(z(1)).^r.*z(2).^(-z(1)*s));

%Run
hw = waitbar(0,'Sampling Gamma Hyper Parameters');
for t = 1:n-1

    %Propose Move
    while(1)
        temp = X(:,t)+SIGMA*randn(2,1);
        if all(temp>=0)
            Y = temp;
            break
        end
    end
    
    %Accept/Reject Move
    ratio = P(Y)*Q(X(:,t),Y,SIGMA)/(P(X(:,t))*Q(Y,X(:,t),SIGMA));
    A = min(1,ratio);
    U = rand;
    if U<=A
        X(:,t+1) = Y;
    else
        X(:,t+1) = X(:,t);
        reject(t+1) = true;
    end
    
    %Update Waitbar
    waitbar(t/n,hw,'Sampling Gamma Hyper Parameters');
     
end
delete(hw)

%Configure Output
a = X(1,:);
b = X(2,:);

