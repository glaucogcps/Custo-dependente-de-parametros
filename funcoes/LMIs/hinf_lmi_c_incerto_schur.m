function out = hinf_lmi_c_incerto_schur(A, B, C, D, varargin)
% function out = hinf_lmi_c_incerto_schur(A, B, C, D, varargin)
%
% Avalia a norma H infinito de um SISTEMA LINEAR CONTÍNUO
% com incerteza politópica, utilizando o Bounded Real Lemma (Schur).
%
% input:
%   A, B, C, D -> Cell arrays (vértices do politopo).
%   varargin   -> Opções: 'solver', 'grau' (de P), 'grau_mu' (do custo).
%
% output:
%   out.hinf   -> Norma Hinf garantida.
%   out.feas   -> 1 se viável, 0 se inviável.
%
% Date: 29/01/2026

%   Bloco de Validação de Entrada  
if nargin &lt; 4
    error('hinf_lmi_c:NotEnoughInputs', 'São necessários A, B, C, D.');
end

input_cells = {A, B, C, D};
num_vertices = numel(A); 

% Validações básicas
n = size(A{1}, 1); m = size(B{1}, 2); p = size(C{1}, 1);
dims_A = size(A{1}); dims_B = size(B{1}); dims_C = size(C{1}); dims_D = size(D{1});

for k = 1:num_vertices
    if ~isequal(size(A{k}), dims_A) || ~isequal(size(B{k}), dims_B) || ...
       ~isequal(size(C{k}), dims_C) || ~isequal(size(D{k}), dims_D)
        error('hinf_lmi_c:VaryingDims', 'Dimensões variam entre vértices.');
    end
end

%   Tratamento das Opções  
options = [];
if nargin > 4
    if nargin == 5 && isstruct(varargin{1})
        options = varargin{1};
    else
        try options = struct(varargin{:}); catch, end
    end
end

if ~isfield(options, 'solver'), options.solver = 'mosek'; end
if ~isfield(options, 'grau'), options.grau = 1; end
if ~isfield(options, 'grau_mu'), options.grau_mu = 0; end
if ~isfield(options, 'verbose'), options.verbose = 0; end

%   Início da Lógica  
out.cpusec_m = clock;

% Variáveis ROLMIP (Sistema)
A_rol = rolmipvar(A, 'A', num_vertices, 1);
B_rol = rolmipvar(B, 'B', num_vertices, 1);
C_rol = rolmipvar(C, 'C', num_vertices, 1);
D_rol = rolmipvar(D, 'D', num_vertices, 1);

% Variável de Lyapunov P(alpha)
P = rolmipvar(n, n, 'P', 'symmetric', num_vertices, options.grau);

% Variável de Custo (gamma^2)
if options.grau_mu == 0
    % Custo constante (gamma^2 escalar)
    gamma2 = sdpvar(1, 1);
    obj = gamma2;
    
    % Termo -gamma^2 * I
    GammaTerm = -gamma2 * eye(m);
    Gamma_rol = rolmipvar(GammaTerm, 'Gamma', num_vertices, 0);
else
    % Custo dependente de parâmetros (mu(alpha))
    % Otimiza o pico de mu(alpha)
    rho = sdpvar(1, 1); % Variável auxiliar para o pico
    
    % mu(alpha) é um polinômio
    mu = rolmipvar(1, 1, 'mu', 'full', num_vertices, options.grau_mu);
    
    % Restrição: mu(alpha) &lt;= rho (para todo alpha)
    % Implementada via coeficientes (condição suficiente)
    coeffs_mu = coeffs(mu);
    LMI_rho = [];
    for i = 1:length(coeffs_mu)
        LMI_rho = [LMI_rho, coeffs_mu{i} &lt;= rho];
    end
    
    obj = rho;
    
    % Termo -mu(alpha) * I
    % Nota: rolmipvar suporta multiplicação escalar-matriz
    Gamma_rol = -mu * eye(m);
end

%   Montagem da LMI (Bounded Real Lemma Contínuo)  
% [ A'P + PA    PB + C'D    C'   ]
% [ B'P + D'C   D'D - g2*I  D'   ] &lt; 0  (Forma incorreta comum)
%
% Forma Correta (Schur Complement do termo D'D):
% [ A'P + PA    PB    C' ]
% [ B'P        -g2*I  D' ] &lt; 0
% [ C           D    -I  ]

% Bloco 1,1: A'P + PA
M11 = A_rol' * P + P * A_rol;

% Bloco 1,2: PB
M12 = P * B_rol;

% Bloco 1,3: C'
M13 = C_rol';

% Bloco 2,2: -gamma^2 * I (ou -mu(alpha)*I)
M22 = Gamma_rol;

% Bloco 2,3: D'
M23 = D_rol';

% Bloco 3,3: -I
M33 = -eye(p);

% Montagem da LMI Simétrica
% [ M11  M12  M13 ]
% [ *    M22  M23 ] &lt; 0
% [ *    *    M33 ]

LMI_sys = [M11, M12, M13;
           M12', M22, M23;
           M13', M23', M33];

LMIs = [P >= 0, LMI_sys &lt;= 0];

if options.grau_mu > 0
    LMIs = [LMIs, LMI_rho];
end

% Informações e Solução
info = lmiinfo(LMIs);
out.L = sum(info.sdp(:, 1));
out.V = size(getvariables(LMIs), 2);
out.cpusec_m = etime(clock, out.cpusec_m);

sol = solvesdp(LMIs, obj, sdpsettings('verbose', options.verbose, 'solver', options.solver));
out.cpusec = sol.solvertime;

% Verificação de viabilidade
p_val = min(checkset(LMIs));
out.p = p_val;

out.hinf = 0;
out.feas = 0;
out.P = [];

if p_val > -1e-6
    out.feas = 1;
    if options.grau_mu == 0
        out.hinf = sqrt(double(gamma2));
    else
        % Se mu é polinomial, o custo garantido é a raiz do pico (rho)
        out.hinf = sqrt(double(rho));
        % Opcional: retornar o objeto mu para análise posterior
        out.mu = mu; 
    end
    out.P = double(P);
end

end