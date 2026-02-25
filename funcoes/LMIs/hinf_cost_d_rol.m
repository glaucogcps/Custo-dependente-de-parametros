function output = hinf_cost_d_rol(A,B,C,D,deg,degGamma)


N=size(A,2);
n=size(A{1},1);
m=size(B{1},2);
p=size(C{1},1);

output.cpusec_m = clock;
Ao = A;
Bo = B;
Co = C;
Do = D;

A = rolmipvar(A, 'A', N, 1);
B = rolmipvar(B, 'B', N, 1);
C = rolmipvar(C, 'C', N, 1);
D = rolmipvar(D, 'D', N, 1);

%generate the variables
P = rolmipvar(n,n, 'P', 'symmetric', N, deg);
op=0; % se op=0 rho tem grau, se op=1 rho é constante
if op == 1
    rho = sdpvar();
    rho_alpha = rolmipvar(1,1, 'rho', 'symmetric', N, degGamma);
    obj = rho;
else        
    powers = gen_coefs(N,degGamma);
    obj = [];
    for i=1:size(powers,1)
        rhos{i} = sdpvar();
        c=factorial(degGamma)/prod(factorial(powers(i,:)));
        if isempty(obj)
            obj = rhos{i}*inv(c);
        else
            obj = rhos{i}*inv(c) + obj;
        end
    end    
    rho_alpha = rolmipvar(rhos, 'rho', N, degGamma);
    %obj = rho_05+0.2*rho_14+0.1*rho_23+0.1*rho_32+0.2*rho_41+rho_50;
    %obj = rhos{1}+0.2*rhos{2}+0.1*rhos{3}+0.1*rhos{4}+0.2*rhos{5}+rhos{6};
    %obj = rho_05+rho_50;
end

varFolga = 1;
if varFolga == 0
    T11 = A'*P+P*A+C'*C;
    T21 = B'*P+D'*C;
    T22 = D'*D-rho_alpha*eye(m);
    T = [T11 T21';T21 T22];
    
    LMIs = [T<=0,P>=0];
    %LMIs = [T<=0,P>=0,[1 rho_alpha;rho_alpha 1]<=rho*eye(2)];
else
    X = rolmipvar(2*n+p+m,n+p, 'X', 'full', N, deg);
    Z = zeros(n);
    Q = blkdiag([P Z;Z -P],eye(p),-rho_alpha*eye(m));
    Bcal = [eye(n) -A zeros(n,p) -B;
            zeros(p,n) -C eye(p) -D];
    LMIs = [Q+X*Bcal+Bcal'*X'<=0,P>=0];
end

% determine the number of LMI rows
info=lmiinfo(LMIs);
% determine the number of LMI rows
output.L =sum(info.sdp(:,1));
% determine the number of scalar variables
output.V=size(getvariables(LMIs),2);

% evaluate the elapsed time to mount the LMIs set
output.cpusec_m = etime(clock,output.cpusec_m);

% solve the LMIs
sol = solvesdp(LMIs,obj,sdpsettings('verbose',0,'solver','sedumi'));

% evaluate the elapsed time to solve the LMIs set
output.cpusec_s =  sol.solvertime;

% retrieving the minimal primal residual
delta=min(checkset(LMIs));

output.feas = 0;
% capturing the solutions (if ones exist)
if delta> -1e-5
    output.feas = 1;
    %output.hinf = sqrt(double(rho));
    maxNorm=-1e10;
    output.alpha = [];
    output.gcosts = [];
    output.realCosts = [];
    for i=1:N
        for j=i+1:N
            for a=0:0.01:1
                s = zeros(1,N);
                s(1,i) = a;
                s(1,j) = 1-a;
                output.alpha = [output.alpha; a];
                v = double(sqrt(evalpar(rho_alpha,{s})));
                output.gcosts = [output.gcosts;v];
                if v > maxNorm
                    maxNorm = v;
                    wc = zeros(1,N);
                    wc(1,i)= a;
                    wc(1,j) = 1-a;
                end
                Aa = a*Ao{i}+(1-a)*Ao{j};
                Ba = a*Bo{i}+(1-a)*Bo{j};
                Ca = a*Co{i}+(1-a)*Co{j};
                Da = a*Do{i}+(1-a)*Do{j};
                s=ss(Aa,Ba,Ca,Da,-1);
                output.realCosts = [output.realCosts;norm(s,inf)];                
            end
        end
    end
    output.wc_alpha=wc;
    output.wc = maxNorm; 
    Aa = zeros(n);
    Ba = zeros(n,m);
    Ca = zeros(p,n);
    Da=zeros(p,m);
    for i=1:N
        Aa = Aa + wc(1,i)*Ao{i};
        Ba = Ba + wc(1,i)*Bo{i};
        Ca = Ca + wc(1,i)*Co{i};
        Da = Da + wc(1,i)*Do{i};
    end
    s=ss(Aa,Ba,Ca,Da);
    output.normSysWc = norm(s,inf);
end

