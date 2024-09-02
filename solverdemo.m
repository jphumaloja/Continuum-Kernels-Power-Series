% This script demonstrates the kernelsolver routine

% define symbolic variables and parameter functions
syms x z y h % the same symbolic variables are used in kernelsolver
Lf = 1; % lambda: function of x and y
Mf = 1; % mu: function of x
Sf = x.*(x+1).*(y-0.5).*x.^2.*(h-0.5); % sigma: function of x, y, and h
Wf = x.*(x+1).*(y-0.5).*exp(x); % W: function of x and y
THf = -70*exp(x*35/pi^2).*y.*(y-1); % theta: function of x and y
Qf = cos(2*pi*y); % q: function of y

% choose parameters for the kernelsolver routine
N = 15; % approximation order (positive integer)
Ny = 2; % reduced order in y (positive integer, not larger than N)
msg = true; % display messages (true/false)
qe = true; % use exact q (true/false)
cfs = false; % seek for exact solution (true/false)
% solve kernels
[K,KB,rd,ef] = kernelsolver(Lf,Mf,Sf,Wf,THf,Qf,N,Ny,msg,qe,cfs);
% Exact solution can be found for these parameter functions, so the exact 
% solution option (cfs) has been set to false to test the power series 
% method. In order to test the exact solution search, set the cfs option to 
% true.
%% plot the control gains
figure(1)
ezsurf(subs(K,x,1),[0 1]);
shading interp
title('')
figure(2)
ezplot(subs(KB,x,1),[0 1])
title('')
% MATLAB complains about ezsurf and ezplot, but they are much faster than
% the recommended fsurf and fplot (or fimplicit)