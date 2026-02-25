function output = h2_lmi_finsler_c_classico(A,B,C,varargin)
% function output = h2_lmi_finsler_c(A,B,C,varargin)
%
% Evaluate the H-2 norm of a continuous-time linear system using extended conditions by Finsler's Lemma. 
% The LMIs are programmed using YALMIP and can be solved by any LMI solver supported by YALMIP (SeDuMi is 
% the default). 
%
% input:  (A,B,C) -> state-space matrices
%        
% output: output.h2       -> H-2 norm (0 if unfeasible)
%         output.cpusec   -> cpu time to solve the LMIs (seconds)
%         output.cpusec_m -> cpu time to mount the LMIs (seconds)
%         output.P        -> Lyapunov matrix
%	      output.V        -> number of decision variables
%	      output.L        -> number of LMIs rows
%
% Date: 12/09/2011
% Author: ricfow@dt.fee.unicamp.br
%
% A = [-0.9  0.2;
%       -0.5  -1.9 ];
% B = [1;0];
% C = [1 0];
% output = h2_lmi_finsler_c(A,B,C)

if nargin > 3
    if nargin == 4
        options = varargin{1};
    else
        options = struct(varargin{:});
    end
else
    options = [];
end
if ~isfield(options,'solver')
    options.solver = 'sedumi';
end

if ~isfield(options,'h2')
    options.h2 = 0;
end
if ~isfield(options,'tolerance')
    options.tolerance = 1e-7;
end
if ~isfield(options,'matLyapDefPos')
    options.matLyapDefPos = 1;
end

%determine the number of vertices of polytope
order = size(A,1);
inputs = size(B,2);
outputs = size(C,1);

output.cpusec_m = clock;

%new LMI system
LMIs = [];

%create the variables
if options.h2 == 0
    mu = sdpvar(1);
    obj = mu;
else
    mu=options.h2*options.h2;
    options.tolerance = 0;
    obj = [];
end
X = sdpvar(outputs,outputs,'symmetric');
H = sdpvar(outputs,order,'full');
J = sdpvar(order,order,'full');
F = sdpvar(order,order,'full');
G = sdpvar(order,order,'full');
% Lyapunov matrix
P = sdpvar(order,order,'symmetric');
if options.matLyapDefPos == 1
    LMIs = [LMIs, P >= 0];
end

%trace condition (mu > trace( B' P B) )
LMIs = [LMIs, mu >= trace(X)];

T11 = C*H'+H*C'-X;
T12 = C*J-H;
T22 = P-J-J';
LMIs = [LMIs, [T11 T12;T12' T22] <= 0];

%gramian condition (A P + P A' + C'C < 0)
T11 = A*F'+F*A';
T12 = P + A*G-F;
T13 = B;
T22 = -G-G';
T23 = zeros(order,inputs);
T33 = -eye(inputs);

LMIs = [LMIs, [T11 T12 T13;T12' T22 T23;T13' T23' T33] <= 0];

output.cpusec_m = etime(clock,output.cpusec_m);

% determine the number of LMI rows
output.L = 0;
for i=1:size(LMIs,1)
    output.L = output.L + size(sdpvar(LMIs(i)),1);
end
% determine the number of scalar variables
output.V = size(getvariables(LMIs),2);

sol = solvesdp(LMIs,obj,sdpsettings('verbose',0,'solver',options.solver));
output.cpusec = sol.solvertime;
p=min(checkset(LMIs));
output.delta = p;

output.feas = 0;
%capturing the solutions (if ones exist)
if p > -options.tolerance
    output.feas = 1;
    output.h2 = sqrt(double(mu));    
    output.P = double(P);
end