function out = h2_lmi_c_incerto_param(A, B, C, D, varargin)
% function out = h2_lmi_c_incerto_param(A, B, C, D, varargin)
%
% Avalia o custo H2 garantido de um SISTEMA LINEAR CONTÍNUO.
% Versão Corrigida e Numericamente Robusta.
%
% Date: 17/02/2026

%   Validação  
if nargin < 4, error('Faltam argumentos'); end
num_vertices = numel(A); 
n = size(A{1}, 1); m = size(B{1}, 2); p = size(C{1}, 1);

%   Opções  
options = struct('solver', 'mosek', 'deg', 1, 'degGamma', 0, 'op', 0, 'varFolga', 1, 'verbose', 0);
if nargin > 4 && isstruct(varargin{1})
    f = fieldnames(varargin{1});
    for i = 1:length(f), options.(f{i}) = varargin{1}.(f{i}); end
end

%   Variáveis ROLMIP  
A_rol = rolmipvar(A, 'A', num_vertices, 1);
B_rol = rolmipvar(B, 'B', num_vertices, 1);
C_rol = rolmipvar(C, 'C', num_vertices, 1);

P = rolmipvar(n, n, 'P', 'symmetric', num_vertices, options.deg);

%   Variável de Custo W(alpha)  
if options.op == 1
    W = rolmipvar(m, m, 'W', 'symmetric', num_vertices, 0); 
    obj = trace(W);
else
    W = rolmipvar(m, m, 'W', 'symmetric', num_vertices, options.degGamma);
    obj = trace(W); 
end

%   LMIs  
if options.varFolga == 0
    %   Clássica (Bounded Real Lemma H2)  
    % 1. Estabilidade: A'P + PA + C'C < 0
    M11 = A_rol' * P + P * A_rol;
    M12 = C_rol';
    M22 = -eye(p);
    LMI_stab = [M11, M12; M12', M22] <= -1e-6*eye(n+p); % Estrita
    
    % 2. Custo: B'PB < W  => [W B'P; PB P] > 0
    N11 = W;
    N12 = B_rol' * P;
    N22 = P;
    LMI_cost = [N11, N12; N12', N22] >= 1e-6*eye(m+n); % Estrita
    
    LMIs = [LMI_stab, LMI_cost];
else
    %   Relaxada (Finsler / OP08b)  
    % Variável de folga G (nxn)
    G = rolmipvar(n, n, 'G', 'full', num_vertices, options.deg);
    
    % 1. Estabilidade Relaxada (Forma Dual - Mais comum para análise)
    % [ A G + G' A'    P - G + G' A'    B ]
    % [ P - G' + A G   -G - G'          0 ] < 0
    % [ B'             0                -I]
    % *Nota: Para H2, a entrada é B (impulso) e saída é C.
    % A condição de estabilidade é sobre A e C.
    
    % Vamos usar a forma PRIMAL que é mais robusta para análise de estabilidade:
    % [ G'A + A'G    G'B        C'       P - G' + A'G ]
    % [ B'G          -I         D'       B'G          ] ... (D=0)
    %
    % Simplificando para H2 (A, C):
    % [ A'G + G'A    P - G + G'A    C' ]
    % [ P - G' + A'G  -G - G'       0  ] < 0
    % [ C             0             -I ]
    
    R11 = A_rol' * G + G' * A_rol;
    R12 = P - G + A_rol' * G'; % Termo de acoplamento
    R13 = C_rol';
    
    R22 = -G - G';
    R23 = zeros(n, p);
    R33 = -eye(p);
    
    LMI_stab = [R11, R12, R13;
                R12', R22, R23;
                R13', R23', R33] <= -1e-6*eye(2*n+p);
                
    % 2. Custo Relaxado
    % Trace(B' P B) < Trace(W)
    % [ W      B'G'     ]
    % [ GB   G + G' - P ] > 0
    
    S11 = W;
    S12 = B_rol' * G';
    S22 = G + G' - P;
    
    LMI_cost = [S11, S12; S12', S22] >= 1e-6*eye(m+n);
    
    LMIs = [LMI_stab, LMI_cost];
end

%   Solução  
info = lmiinfo(LMIs);
out.L = sum(info.sdp(:, 1));
out.V = size(getvariables(LMIs), 2);

out.cpusec_m = clock;
sol = solvesdp(LMIs, obj, sdpsettings('verbose', options.verbose, 'solver', options.solver));
out.cpusec_s = sol.solvertime;

if min(checkset(LMIs)) > -1e-6
    out.feas = 1;
    
    % Calcula pior caso (Grid)
    max_val = -inf;
    out.gcosts = [];
    out.realCosts = [];
    out.alpha = 0:0.05:1; % Salva o vetor alpha para plotagem
    
    % Preserva originais para validação
    Ao = A; Bo = B; Co = C;
    
    for a = out.alpha
        % 1. Custo Garantido
        W_val = double(evalpar(W, {[a, 1-a]}));
        v = sqrt(trace(W_val));
        out.gcosts = [out.gcosts; v];
        if v > max_val, max_val = v; end
        
        % 2. Custo Real
        Aa = a*Ao{1} + (1-a)*Ao{2};
        Ba = a*Bo{1} + (1-a)*Bo{2};
        Ca = a*Co{1} + (1-a)*Co{2};
        sys = ss(Aa, Ba, Ca, zeros(p,m), 0);
        try
            out.realCosts = [out.realCosts; norm(sys, 2)];
        catch
            out.realCosts = [out.realCosts; NaN];
        end
    end
    out.wc = max_val;
else
    out.feas = 0;
    out.wc = Inf;
    out.gcosts = [];
    out.realCosts = [];
    out.alpha = [];
end
end