function [Ks,KBs,rd,ef] = kernelsolver(Lf,Mf,Sf,Wf,THf,Qf,N,Ny,msg,qe,cfs)
%KERNELSOLVER Solves continuum kernel equations for given parameters
%   Output arguments:
%   Ks: ensemble part of the continuum solution
%   KBs: +1 part of the continuum solution
%   rd: residual of the least squares fit of the K, KB coefficients
%   ef: exit flag (=0 for exact, >0 for power series solution)
%   Input arguments
%   Lf: lambda (symbolic function of x and y)
%   Mf: mu (symbolic function of x)
%   Sf: sigma (symbolic function of x, y and h)
%   Wf: W (symbolic function of x and y)
%   THf: theta (symbolic function of x and y)
%   Qf: q (symbolic function of y)
%   N: order of the power series approximation
%   Ny: reduced approximation order in y (optional, default: N)
%   msg: true/false to display messages (optional, default: true)
%   qe: true/false to use exact q (optional, default: false)
%   cls: true/false to seek exact solution (optional, default: false)

% use the same symbolic variables as in the script
syms x z y h
% checks for input arguments
if nargin < 7
  error('Too few input arguments')
end
if nargin < 8
  Ny = N;
end
if nargin < 9
  msg = true;
end
if nargin < 10
  qe = false;
end
if nargin < 11
  cfs = false;
end
if N <= 0
  error('The order N must be positive')
end
if N < Ny
  error('The order Ny cannot be larger than N')
end
% seek for exact solution (if selected)
if cfs
  [Ks,KBs,ef] = closedform(Lf,Mf,Sf,Wf,THf,Qf,msg);
  if ~ef % if exact solution found, set residual to 0 and return
    if msg
      disp('Exact solution found')
    end
    rd = 0;
    return
  else
    if msg
      disp('Could not find exact solution, initiating power series')
    end
  end
end

% compute Taylor series for parameters, except use exact q if selected
L = taylor(Lf,[x y],'order',N+1);
M = taylor(Mf,x,'order',N+1);
W = taylor(Wf,[x y],'order',N+1);
S = taylor(Sf,[x y h],'order',N+1);
TH = taylor(THf,[x y],'order',N+1);
if qe
  Q = Qf;
else
  Q = taylor(Qf,y,'order',N+1);
end
% initialize power series for K and KB
NK = N*(N+1)*(2*N+10)/12+N+1; % number of terms for K (full)
if Ny < N % if reduced order in y used
  NY = N-Ny-1; % aux variable
  NKy = NY*(NY+1)*(2*NY+10)/12+NY+1; % number of reduced terms
  NK = NK-NKy; % numer of terms for K (reduced)
end
NB = (N+1)*(N+2)/2; % number of terms for Kbar
KC = sym('k', [1 NK+NB],'real'); % symbols for coefficients
IK = repmat("0",1,NK+NB); % tags for coefficients (not used)
rind = 0;
K = 0;
for tk=0:Ny
  for tm=0:N-tk
    for tn=0:tm
      rind = rind + 1;
      K = K + KC(rind)*x^tn*z^(tm-tn)*y^tk;
      IK(rind) = [num2str(tn),num2str(tm-tn),num2str(tk)];
    end
  end
end
% Kbar Taylor (2D)
KB = 0;
for tm=0:N
  for tn=0:tm
    rind = rind + 1;
    KB = KB + KC(rind)*x^tn*z^(tm-tn);
    IK(rind) = [num2str(tn) num2str(tm-tn)];
  end
end
% construct equations and boundary conditions with power series
KE1 = M*diff(K,x) -  subs(L,x,z)*diff(K,z) ... 
  - subs(TH,x,z)*KB - subs(diff(L,x),x,z)*K ...
  + int(subs(S,[x,y,h],[z,h,y])*subs(K,y,h),h,0,1);
KE2 = M*diff(KB,x) + subs(M,x,z)*diff(KB,z) ...
  + subs(diff(M,x),x,z)*KB - int(subs(W,x,z)*K,y,0,1);
BC1 = (L+M)*subs(K,z,x) + TH;
BC2 = subs(M,x,0)*subs(KB,z,0) - int(Q*subs(L,x,0)*subs(K,z,0),y,0,1);
% extract coefficients for different powers of spatial variables
KE1C = coeffs(KE1,[x z y]);
KE2C = coeffs(KE2,[x z]);
BC1C = coeffs(BC1,[x y]);
BC2C = coeffs(BC2,x);
% construct set linear equations for the coefficients based on the above%
neqs = [numel(KE1C), numel(KE2C), numel(BC1C), numel(BC2C)];
A = zeros(sum(neqs), NK+NB); % initialize A...
B = zeros(sum(neqs),1); % ...and b
EQS = {KE1C,KE2C,BC1C,BC2C}; % cell array for the equation data
if msg
  disp(['Parsing data: ',num2str(sum(neqs)),' equations for ',...
    num2str(NK+NB),' unknonws']);
end
rind = 0; % row index
for ii = 1:4
  DAT = EQS{ii};
  for m=1:neqs(ii)
    rind = rind + 1;
    [cc, kc] = coeffs(DAT(m),KC);
    for mk=1:numel(kc) % check which K's are present in the coefficients
      str = char(kc(mk)); 
      sl = numel(str);
      if sl > 1 % look up the index and insert to A
        A(rind,str2double(str(2:sl))) = cc(mk);
      else % if no index, it's a contant; insert to B
        B(rind) = -cc(mk);
      end
    end
  end
  if msg
    switch ii
      case 1
        disp('First kernel equation parsed')
      case 2
        disp('Second kernel equation parsed')
      case 3
        disp('First boundary condtition parsed')
      case 4
        disp('Second boundary condition parsed')
    end
  end
end
% solve equations and compute residual
Kc = A\B;
rd = norm(A*Kc-B);
if ~cfs
  ef = 8;
end
% insert the obtained values into the power series
Ks = subs(K,KC(1:NK),Kc(1:NK)');
KBs = subs(KB,KC(NK+1:NK+NB),Kc(NK+1:NK+NB)');
end

function [Ks, KBs, ef] = closedform(Lf,Mf,Sf,Wf,THf,Qf,msg)
%CLOSEDFORM Seeks exact solution to continuum kernel equations
%   Parameters are the same as for KERNELSOLVER
%   ef = 0 for success, ef > 0 for failure (check messages for reason)

% use the same symbolic variables as in the main function and script
syms x z y h
% if any check fails, return empty K and KB, and the corresponding ef
Ks = [];
KBs = [];
% test conditions for explicit solution
if has(sym(Lf),[x y])
  if msg
    disp('The parameter lambda is not constant')
  end
  ef = 1;
  return
end
if has(sym(Mf),x)
  if msg
      disp('The parameter mu is not constant')
  end
  ef = 2;
  return
end
% check separability of sigma, W, and theta, starting with sigma
% find some integers for which sigma is not zero
for kk = 1:100
  vec = randi(5,1,3);
  if ~isequal(expand(subs(Sf,[x y h],vec)),sym(0))
    break
  end
  if kk == 100 % if not found, return an error
    error('Could not check if sigma is separable')
  end
end
% test if sigma is separable (if yes, these are the components)
Sfx = subs(Sf,[y h], vec(2:3))/subs(Sf,[x y h],vec);
Sfy = subs(Sf,[x h], vec([1 3]))/subs(Sf,[x y h],vec);
Sfh = subs(Sf,[x y], vec(1:2));
if ~isequal(expand(Sfx*Sfy*Sfh),expand(Sf))
  if msg
    disp('The parameter sigma is not separable')
  end
  ef = 3;
  return
end
% find some integers for which W is not zero
for kk = 1:100
  vec = randi(5,1,2);
  if ~isequal(expand(subs(Wf,[x y],vec)),sym(0))
    break
  end
  if kk == 100 % if not found, return an error
    error('Could not check if W is separable')
  end
end
% test if W is separable (if yes, these are the components)
Wfx  = subs(Wf,y,vec(2))/subs(Wf,[x y],vec);
Wfy = subs(Wf,x,vec(1));
if ~isequal(expand(Wfx*Wfy),expand(Wf))
  if msg
    disp('The parameter W is not separable')
  end
  ef = 4;
  return
end
% find some integers for which TH is not zero
for kk = 1:100
  vec = randi(5,1,2);
  if ~isequal(expand(subs(THf,[x y],vec)),sym(0))
    break
  end
  if kk == 100 % if not found, return an error
    error('Could not check if theta is separable')
  end
end
% test if TH is separable (if yes, these are the components)
THfx  = subs(THf,y,vec(2))/subs(THf,[x y],vec);
THfy = subs(THf,x,vec(1));
if ~isequal(expand(THfx*THfy),expand(THf))
  if msg
    disp('The parameter theta is not separable')
  end
  ef = 5;
  return
end
% compute c_y, test if it is a constant
cy = simplify(subs(Sfh,h,y)/THfy*int(Sfy*THfy,[0 1]));
if has(cy,y)
  if msg
    disp('the constant c_y does not exist')
  end
  ef = 6;
  return
end
% check condition (7)
lhs = cy*diff(Sfx) + simplify(Lf*(diff(THfx,2)*THfx-diff(THfx)^2)/THfx^2);
rhs = Wfx*THfx*int(Wfy*THfy,y,[0 1]);
if ~isequal(expand(lhs),expand(rhs))
  if msg
    disp('the condition (7) is not satisfied')
  end
  ef = 7;
  return
end
% compute the explicit solution based on (8)-(10)
cx = Mf/(Lf+Mf)*(cy*subs(Sfx,x,0) + ...
  Lf*simplify(subs(diff(THfx),x,0)/subs(THfx,x,0)) +...
  Lf/Mf*subs(THfx,x,0)*int(Qf*THfy,[0 1]));
KBz = (Lf*(simplify(subs(diff(THfx),x,z)/subs(THfx,x,z)) - cx/Mf) + ...
  cy*subs(Sfx,x,z) - cx)*exp(-cx/Mf*z)/(Lf+Mf);
Ks = -1/(Lf+Mf)*exp(cx/Mf*x)*exp(-cx/Mf*z)*subs(THf,x,z);
KBs = exp(cx/Mf*x)*KBz;
ef = 0; % exit flag 0 for exact solution
end