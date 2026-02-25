function out = hinf_lmi_c_incerto_param(A, B, C, D, varargin)
% function out = hinf_lmi_c_incerto_param(A, B, C, D, varargin)
%
% Avalia o custo H infinito garantido de um SISTEMA LINEAR CONTÍNUO.
% Retorna o polinômio de custo para plotagem e estatísticas de complexidade.
%
% Date: 17/02/2026

%   Validação de Entrada  
if nargin < 4, error('Faltam argumentos A, B, C, D'); end
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
D_rol = rolmipvar(D, 'D', num_vertices, 1);

P = rolmipvar(n, n, 'P', 'symmetric', num_vertices, options.deg);

%   Variável de Desempenho (Gamma^2)  
if options.op == 1
    rho = sdpvar(1, 1);
    rho_alpha = rolmipvar(rho, 'rho', num_vertices, 0);
    obj = rho;
else        
    powers = gen_coefs(num_vertices, options.degGamma);
    rhos = cell(1, size(powers, 1));
    obj = 0;
    for i = 1:size(powers, 1)
        rhos{i} = sdpvar(1, 1);
        c_val = factorial(options.degGamma) / prod(factorial(powers(i, :)));
        obj = obj + rhos{i} * (1/c_val);
    end    
    rho_alpha = rolmipvar(rhos, 'rho', num_vertices, options.degGamma);
end

%   LMIs  
if options.varFolga == 0
    % Bounded Real Lemma Clássico
    M11 = A_rol' * P + P * A_rol;
    M12 = P * B_rol;
    M13 = C_rol';
    M22 = -rho_alpha * eye(m);
    M23 = D_rol';
    M33 = -eye(p);
    
    LMI_sys = [M11, M12, M13; M12', M22, M23; M13', M23', M33];
    LMIs = [LMI_sys <= 0, P >= 1e-6*eye(n)];
else
    % Finsler / Relaxada
    X = rolmipvar(2*n + p + m, n + p, 'X', 'full', num_vertices, options.deg);
    Z_n = zeros(n);
    
    Q = blkdiag([Z_n, P; P, Z_n], eye(p), -rho_alpha * eye(m));
    
    Bcal = [eye(n), -A_rol, zeros(n, p), -B_rol;
            zeros(p, n), -C_rol, eye(p), -D_rol];
            
    LMIs = [Q + X * Bcal + Bcal' * X' <= 0, P >= 1e-6*eye(n)];
end

%   Estatísticas de Complexidade (CORREÇÃO AQUI)  
info = lmiinfo(LMIs);
out.L = sum(info.sdp(:, 1));          % Número total de linhas de LMI
out.V = size(getvariables(LMIs), 2);  % Número total de variáveis escalares

%   Solução  
out.cpusec_m = clock;
sol = solvesdp(LMIs, obj, sdpsettings('verbose', options.verbose, 'solver', options.solver));
out.cpusec_s = sol.solvertime;

if min(checkset(LMIs)) > -1e-6
    out.feas = 1;
    out.rho_alpha = rho_alpha; 
    
    % Calcula o pior caso (máximo do polinômio no simplex)
    max_val = -inf;
    for a = 0:0.05:1
        val = double(evalpar(rho_alpha, {[a, 1-a]}));
        if val > max_val, max_val = val; end
    end
    out.wc = sqrt(max(0, max_val)); 
else
    out.feas = 0;
    out.rho_alpha = [];
    out.wc = Inf;
end
end