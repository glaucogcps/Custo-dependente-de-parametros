function out = hinf_lmi_c_incerto_schur(A, B, C, D, varargin)
% function out = hinf_lmi_c_incerto_schur(A, B, C, D, varargin)
%
% Avalia o custo H infinito garantido de um SISTEMA LINEAR CONTÍNUO
% com incerteza politópica, utilizando o Bounded Real Lemma (Schur).
%
% Entradas:
%   A, B, C, D - Cell arrays contendo as matrizes dos vértices do politopo.
%                Ex: A = {A1, A2, ..., AN}, onde N é o número de vértices.
%   varargin   - (Opcional) Struct ou pares nome-valor com opções:
%                .solver   : Solver a ser utilizado (default: 'mosek')
%                .deg      : Grau do polinômio da matriz de Lyapunov P (default: 1)
%                .degGamma : Grau do polinômio da variável de desempenho mu (default: 0)
%                .verbose  : Nível de detalhe do solver (default: 0)
%                *Nota: 'grau' e 'grau_mu' são aceitos por retrocompatibilidade.
%
% Saída:
%   out - Struct contendo os resultados da otimização e validação:
%         .feas      : Status de viabilidade (1 = viável, 0 = inviável)
%         .hinf      : Pior custo Hinf garantido (pico) retornado pelo solver
%         .wc        : Pior custo Hinf garantido encontrado no grid
%         .normSysWc : Custo Hinf real do sistema no ponto de pior caso garantido
%         .wc_alpha  : Coordenadas baricêntricas do pior caso
%         .gcosts    : Vetor com os custos garantidos avaliados no grid
%         .realCosts : Vetor com os custos reais (norma Hinf) avaliados no grid
%         .alpha     : Matriz com os pontos do simplex testados (cada linha é um ponto)
%         .L         : Número de linhas das LMIs
%         .V         : Número de variáveis escalares
%         .cpusec_m  : Tempo de montagem das LMIs
%         .cpusec_s  : Tempo de resolução do solver
%
% Date: 22/03/2026

%   Bloco de Validação de Entrada  
if nargin < 4
    error('hinf_lmi_c:NotEnoughInputs', 'São necessários A, B, C, D (cell arrays de matrizes).');
end

input_cells = {A, B, C, D};
input_names = {'A', 'B', 'C', 'D'};
num_vertices = numel(A); % N

% Valida tipos e consistência do número de vértices
for i = 1:length(input_cells)
    if ~iscell(input_cells{i})
        error('hinf_lmi_c:InvalidInputType', 'O argumento %s deve ser um cell array.', input_names{i});
    end
    if numel(input_cells{i}) ~= num_vertices
        error('hinf_lmi_c:DimensionMismatch', 'Todos os cell arrays devem ter o mesmo número de vértices.');
    end
end

% Valida dimensões das matrizes
n = size(A{1}, 1); m = size(B{1}, 2); p = size(C{1}, 1);
dims_A = size(A{1}); dims_B = size(B{1}); dims_C = size(C{1}); dims_D = size(D{1});

if dims_A(1) ~= dims_A(2)
    error('hinf_lmi_c:InvalidADim', 'Matrizes A devem ser quadradas.');
end
if dims_B(1) ~= n || dims_C(2) ~= n
    error('hinf_lmi_c:IncompatibleDims', 'Dimensões de B ou C incompatíveis com A.');
end
if dims_D(1) ~= p || dims_D(2) ~= m
    error('hinf_lmi_c:IncompatibleDims', 'Dimensões de D incompatíveis com C e B.');
end

% Valida conteúdo numérico
for k = 1:num_vertices
    if ~isnumeric(A{k}) || ~isnumeric(B{k}) || ~isnumeric(C{k}) || ~isnumeric(D{k})
        error('hinf_lmi_c:NonNumeric', 'Conteúdo das células deve ser numérico.');
    end
    if ~isequal(size(A{k}), dims_A) || ~isequal(size(B{k}), dims_B) || ...
       ~isequal(size(C{k}), dims_C) || ~isequal(size(D{k}), dims_D)
        error('hinf_lmi_c:VaryingDims', 'Dimensões das matrizes variam entre os vértices.');
    end
end

%   Tratamento das Opções  
options = struct('solver', 'mosek', 'deg', 1, 'degGamma', 0, 'verbose', 0);
if nargin > 4
    if nargin == 5 && isstruct(varargin{1})
        opts_in = varargin{1};
    else
        try opts_in = struct(varargin{:}); catch, opts_in = struct(); end
    end
    
    f = fieldnames(opts_in);
    for i = 1:length(f)
        % Mapeamento de retrocompatibilidade
        if strcmp(f{i}, 'grau'), options.deg = opts_in.grau;
        elseif strcmp(f{i}, 'grau_mu'), options.degGamma = opts_in.grau_mu;
        else, options.(f{i}) = opts_in.(f{i}); end
    end
end

%   Início da Lógica  
out.cpusec_m = clock;

% Preserva as matrizes originais para a validação final
Ao = A; Bo = B; Co = C; Do = D;

% Variáveis ROLMIP (Sistema)
A_rol = rolmipvar(A, 'A', num_vertices, 1);
B_rol = rolmipvar(B, 'B', num_vertices, 1);
C_rol = rolmipvar(C, 'C', num_vertices, 1);
D_rol = rolmipvar(D, 'D', num_vertices, 1);

% Variável de Lyapunov P(alpha)
P = rolmipvar(n, n, 'P', 'symmetric', num_vertices, options.deg);

% Variável de Custo (gamma^2)
if options.degGamma== 0
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

if options.degGamma> 0
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
out.hinf = Inf;
out.feas = 0;
out.P = [];
out.wc = Inf;
out.normSysWc = Inf;
out.alpha = [];
out.gcosts = [];
out.realCosts = [];
out.wc_alpha = [];

if p_val > -1e-6
    out.feas = 1;
    
    if options.degGamma == 0
        out.hinf = sqrt(double(gamma2));
    else
        % Se mu é polinomial, o custo garantido é a raiz do pico (rho)
        out.hinf = sqrt(double(rho));
        out.mu = mu; 
    end
    
    %   Validação (Grid Search em todo o Simplex)  
    maxNorm = -1e10;
    passo_grid = 0.05; 
    
    % Gera os pontos do simplex unitário usando a nova função
    pontos_alpha = particao_simplex(num_vertices, passo_grid);
    out.alpha = zeros(length(pontos_alpha), num_vertices);
    
    for k = 1:length(pontos_alpha)
        alpha_vec = pontos_alpha{k}; % Vetor coluna
        s = alpha_vec'; % Vetor linha para o ROLMIP
        out.alpha(k, :) = s;
        
        % 1. Custo Garantido (via LMI)
        if options.degGamma == 0
            v = out.hinf; % Custo é constante
        else
            v = double(sqrt(evalpar(mu, {s}))); % Avalia o polinômio mu(alpha)
        end
        out.gcosts = [out.gcosts; v];
        
        % Rastreia o pior caso garantido
        if v > maxNorm
            maxNorm = v;
            wc_vec = s;
        end
        
        % 2. Custo Real (Sistema Congelado)
        Aa_val = zeros(n, n); Ba_val = zeros(n, m);
        Ca_val = zeros(p, n); Da_val = zeros(p, m);
        
        % Combinação convexa generalizada para N vértices
        for i = 1:num_vertices
            Aa_val = Aa_val + alpha_vec(i) * Ao{i};
            Ba_val = Ba_val + alpha_vec(i) * Bo{i};
            Ca_val = Ca_val + alpha_vec(i) * Co{i};
            Da_val = Da_val + alpha_vec(i) * Do{i};
        end
        
        % Cria sistema SS contínuo
        sys_local = ss(Aa_val, Ba_val, Ca_val, Da_val);
        try
            % CRÍTICO: Usar norma inf para H-infinito
            out.realCosts = [out.realCosts; norm(sys_local, inf)];
        catch
            out.realCosts = [out.realCosts; NaN];
        end
    end
    
    out.wc_alpha = wc_vec;
    out.wc = maxNorm; 
    
    % Recalcula o sistema no ponto de pior caso encontrado
    Aa_wc = zeros(n, n); Ba_wc = zeros(n, m); Ca_wc = zeros(p, n); Da_wc = zeros(p, m);
    for i = 1:num_vertices
        Aa_wc = Aa_wc + wc_vec(i) * Ao{i};
        Ba_wc = Ba_wc + wc_vec(i) * Bo{i};
        Ca_wc = Ca_wc + wc_vec(i) * Co{i};
        Da_wc = Da_wc + wc_vec(i) * Do{i};
    end
    
    sys_wc = ss(Aa_wc, Ba_wc, Ca_wc, Da_wc);
    try
        out.normSysWc = norm(sys_wc, inf);
    catch
        out.normSysWc = NaN;
    end
end

end