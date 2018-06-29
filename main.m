%MAIN() Executive script for running metropolisHastings.m 
%   This script runs metropolisHastings.m on simulated data.
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

%Clean Up
clear
close all
clc

%Set Random Number Generator
SEED = 1;
rng(SEED)

%Generate Simulated Observations
atrue = 1;
btrue = 2;
N = 200;
x = gamrnd(atrue,1/btrue,1,N);

%Compute Posterior Parameters
p = prod(x);
q = sum(x);
r = N;
s = N;

%Evaluate Posterior Surface
gx = 0:0.01:2;
gy = 0:0.01:4;
[A,B] = meshgrid(gx,gy);
F = p.^(A-1).*exp(-B*q)./(gamma(A).^r.*B.^(-A*s));
F(isnan(F)) = inf;

%Run Metropolis-Hastings Sampler
nsamp = 5000;
sigma = 0.1;
icond = [0;0];
[a,b,reject] = metropolisHastings(x,nsamp,'sigma',sigma,'icond',icond);

%Plot Markov Chain Trajectory
figure,
subplot(121)
mesh(gx,gy,F)
xlabel('a'),ylabel('b'),title('Posterior Distribution Heat Map')
view(0,90)
subplot(122)
plot(a,b,'k','linewidth',0.5)
xlabel('a'),ylabel('b'),title('Posterior Distribution Samples')
axis([0,2,0,4])
grid on

%Plot Hyperparameter Samples
figure,
subplot(211)
plot(a),title('alpha'),xlabel('Iteration'),grid on
hold on, ax = axis; plot(ax(1:2),atrue*[1,1],'k--')
legend('MCMC Samples','Truth','location','southeast')
subplot(212)
plot(b),title('beta'),xlabel('Iteration'),grid on
hold on, ax = axis; plot(ax(1:2),btrue*[1,1],'k--')
legend('MCMC Samples','Truth','location','southeast')

%Plot Rejection Percentage
figure
plot(filter(ones(1,100)/100,1,reject)*100)
grid on
xlabel('MCMC Iteration')
ylabel('MCMC Rejection %')








