function output = hinf_norm_c_classico(A,B,C,D,varargin)
%function output = hinf_norm_c(A,B,C,D,varargin)
%
% Compute the H-infininty norm of a continuous-time linear system using linear
% matrix inequalities. The LMIs are  programmed using YALMIP and can be solved 
% by any LMI solver supported by YALMIP (SeDuMi was used with the default options).
% inputs:  A,B,C,D   -> vertices of the polytope (cell array)
%          optional 
%          hinf -> if this parameter is different of zero, the routine will
%                        check the feasibility of the LMIs for this value of the hinf
%                        norm.
%          tol -> tolerance of the feasbility of the LMIs. 
%          solver -> SDP solver to be used to solve the LMIs.
%
% outputs: output.hinf     -> H-infinity guaranteed cost (0 if unfeasible)
%          output.cpusec   -> cpu time to solve the LMIs (seconds)
%          output.cpusec_m -> cpu time to mount the LMIs (seconds)
%          output.Pk       -> Solution variables P
%          output.V        -> number of scalar variables used in the optmization problem
%          output.L        -> number of LMI rows used in the optmization problem
%          output.delta    -> minimal primal residual returned by the LMI solver (SeDuMi is the default).
%
% Example: 2 states and 2 vertices
% A = [-0.9  0.2;
%       -0.5  -1.9 ];
% B = [1;0];
% C = [1 0];
% D = 0;
% output = hinf_norm_c(A,B,C,D)
%
% Date: 29/08/2011
% Author: ricfow@dt.fee.unicamp.br

options = [];
if nargin > 4
    if nargin == 5
        options = varargin{1};
    else
        options = struct(varargin{:});
    end
end

if ~isfield(options,'solver')
    options.solver = 'sedumi';
end
if ~isfield(options,'tol')
    options.tol = 1e-7;
end
if ~isfield(options,'hinf')
    options.hinf = 0;
end

%determine the number of vertices of polytope
order = size(A,1);
inputs = size(B,2);
outputs = size(C,1);

output.cpusec_m = clock;

%new LMI system
LMIs = [];

if options.hinf == 0
    mu = sdpvar(1);
    obj = mu;
else
    mu = hinf^2;
    obj = [];
    options.tol=0;
end

P = sdpvar(order,order,'symmetric');
LMIs = [LMIs, P >= 0];

% Bounded real lemma  < 0
T11 = A'*P + P*A;
T12 = P*B;
T13 = C';
T22 = -eye(inputs);
T23 = D';
T33 = -mu*eye(outputs);
T = [T11 T12 T13;
    T12' T22 T23;
    T13' T23' T33];
LMIs = [LMIs, T <= 0];

output.cpusec_m = etime(clock,output.cpusec_m);

% determine the number of LMI rows
info=lmiinfo(LMIs);
output.L =sum(info.sdp(:,1));% determine the number of scalar variables
% determine the number of scalar variables
output.V = size(getvariables(LMIs),2);

% solve the LMIs
sol = solvesdp(LMIs,obj,sdpsettings('verbose',0,'solver',options.solver));
% evaluate the elapsed time to solve the LMIs set
output.cpusec_s = sol.solvertime;

% retrieving the minimal primal residual
p=min(checkset(LMIs));
output.delta = p;

output.feas = 0;
% capturing the solutions (if ones exist)
if p > -options.tol
    output.P = double(P);
    output.hinf = sqrt(double(mu));
    output.feas = 1;
end