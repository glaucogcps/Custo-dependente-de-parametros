function out = h2_lmi_c_classico(A,B,C,varargin)
% function out = h2_lmi_c(A,B,C,varargin)
%
% Evaluate the H-2 norm of a continuous-time linear system. The LMIs are programmed using YALMIP and can be 
% solved by any LMI solver supported by YALMIP (SeDuMi is the default). 
%
% input:  (A,B,C) -> state-space matrices
%        
% output: out.h2       -> H-2 norm (0 if unfeasible)
%         out.cpusec   -> cpu time to solve the LMIs (seconds)
%         out.cpusec_m -> cpu time to mount the LMIs (seconds)
%         out.P        -> Lyapunov matrix
%	      out.V        -> number of decision variables
%	      out.L        -> number of LMIs rows
%
% Date: 04/09/2009
% Author: ricfow@dt.fee.unicamp.br

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

%determine the number of vertices of polytope
order = size(A,1);
inputs = size(B,2);
outputs = size(C,1);

out.cpusec_m = clock;

%create the variables
if options.h2 == 0
    mu = sdpvar(1);
    obj = mu;
else
    mu=options.h2*options.h2;
    obj = [];
    options.tolerance = 0;
end

P = sdpvar(order,order,'symmetric');
LMIs = [P >=0];
    
%trace condition (mu > trace( B' P B) )
LMIs = [LMIs, mu >= trace(B'*P*B)];

%gramian condition (A P + P A' + C'C < 0)
LMIs = [LMIs, A'*P + P*A + C'*C <= 0];

out.cpusec_m = etime(clock,out.cpusec_m);

% determine the number of LMI rows
output.L = 0;
for i=1:size(LMIs,1)
    output.L = output.L + size(sdpvar(LMIs(i)),1);
end
% determine the number of scalar variables
output.V = size(getvariables(LMIs),2);


sol = solvesdp(LMIs,obj,sdpsettings('verbose',0,'solver',options.solver));
out.cpusec = sol.solvertime;
p=min(checkset(LMIs));

out.h2 = 0;
%capturing the solutions (if ones exist)
if p > -options.tolerance
    out.h2 = sqrt(double(mu));    
    out.P = double(P);
end
