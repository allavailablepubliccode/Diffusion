function LAP = Diffusion_DCM_inversion(data, struct)
% input 'data' is the nxt timeseries
% input 'struct' is the nxn structural adjacency matrix 

% calculate prior variance on functional connectivity matrix so that
% 10x as many time points as free parameters are retained
if (size(data,1)^2+1)*10 < size(data,2)
    prvr                            = ones(size(data,1));
else
    [U,S]                           = svd(data);
    Um                              = U(:,1:round(size(data,1)/3));
    Km                              = kron(Um*Um', Um*Um');
    Sigma                           = ones(size(Km,1));
    prvr                            = Km*Sigma*Km';
    prvr                            = diag(prvr);
    prvr                            = reshape(prvr,size(U,1),size(U,2));
    prvr(logical(eye(size(prvr))))  = min(prvr(:))*(1/2);
    [~,index]                       = sort(prvr(:),'descend');
    index                           = index(1:floor(size(data,2)/10)-1);
    prvr                            = zeros(size(prvr));
    prvr(index)                     = 1;
end

Errval                              = 4; % noise assumption
Varval                              = 1; % prior variance

N                                   = size(data,2); % number of timepoints
n                                   = size(data,1); % number of regions

struct(logical(eye(size(struct))))  = 0; % no self-structural connectivity

% prior means (adjust according to known constraints)
P.Adj                               = struct; % structural connectivity
P.C                                 = ones(n,1)/16; % extrinsic connect. 
P.sig                               = 0; % diffusion coefficient, note that
% models in which the estimate of the diffusion coefficient is negative 
% should be excluded prior to any subsequent model averaging
P.A                                 = full(-speye(n)/16); % func. conn.

% prior variance
pC.Adj                              = zeros(n); % if structure known
pC.C                                = zeros(n,1); % adjust if to be varied
pC.sig                              = 1;
pC.A                                = prvr;

DEM.M(1).x.s                        = ones(n,1); % initial states
DEM.M(1).f                          = 'f_fun(x,v,P)';
DEM.M(1).g                          = @(x,v,P) x.s;
DEM.M(1).pE                         = P;
DEM.M(1).pC                         = diag(spm_vec(pC))*Varval;
DEM.M(1).V                          = exp(Errval);
DEM.M(1).W                          = exp(Errval);

DEM.M(2).v                          = 0;
DEM.M(2).V                          = exp(Errval);
DEM.Y                               = data;
DEM.U                               = rand(1,N)/16; % random external inp.

LAP                                 = spm_DEM(DEM);

% Diffusion_DCM_inversion calls f_fun below
% function f = f_fun(x,v,P)
% 
% Adjmat                          = P.Adj;
% Sig                             = P.sig;
% Amat                            = P.A;
% Amain                           = diag(Amat);
% 
% Jac                             = Amat - Sig*Adjmat;
% 
% for ii                          = 1:size(Amat,1)
% Jmain(ii)                       = Amain(ii) + ...
%                                   Sig*(sum(Adjmat(ii,1:ii-1)) + ...
%                                   sum(Adjmat(ii,ii+1:end)));
% end
% 
% Jac(logical(eye(size(Jac))))    = Jmain;
% 
% f                               = Jac*x.s + P.C*v;
