function output = h2_lmi_mjls_d_incerto_controle(A, B, E, C, D, Gamma, mu, varargin)
% function output = h2_lmi_mjls_d_incerto_controle(A, B, E, C, D, Gamma, mu, varargin)
%
% Synthesize an H-2 guaranteed cost controller (u = K_i x) for a discrete-time 
% MJLS with uncertain probability matrix (polytopic).
%
% input:  (A,B,E,C,D) -> state-space matrices (cell arrays for modes)
%         Gamma       -> cell array with vertices of the probability matrix polytope
%         mu          -> initial probability distribution
%
% varargin (Opções):
%   .solver    : Solver LMI a ser utilizado (default: 'sedumi')
%   .deg       : Grau do polinômio da matriz de Lyapunov P (default: 1)
%   .h2        : Limite superior desejado para a norma H-2. Se for 0, a 
%                função minimiza a norma H-2 (default: 0)
%   .tolerance : Tolerância numérica para a verificação de viabilidade (default: 1e-7)
%   .op        : Tipo de otimização (0 = integral/média, 1 = pico) (default: 1)
%   .degrho    : Grau do polinômio da variável de desempenho rho (default: 0)
%   .indep     : Controlador independente de modo (1 = sim, 0 = não) (default: 0)
%        
% Saída (output):
%   .feas      : Status de viabilidade (1 = viável, 0 = inviável)
%   .h2        : Pior custo H-2 garantido (mantido para retrocompatibilidade)
%   .wc        : Pior custo H-2 garantido encontrado no grid
%   .normSysWc : Custo H-2 real do sistema em malha fechada no pior caso
%   .K         : Ganhos do controlador K_i para cada modo (cell array)
%   .Acl...    : Matrizes do sistema em malha fechada (Acl, Ecl, Ccl, Ezcl)
%   .wc_alpha  : Ponto do simplex onde ocorreu o pior caso
%   .P         : Matrizes de Lyapunov P_i avaliadas
%   .cpusec    : Tempo de solução do solver (segundos)
%   .cpusec_m  : Tempo de montagem das LMIs (segundos)
%   .V         : Número de variáveis de decisão escalares
%   .L         : Número de linhas das LMIs
%   .gcosts    : Vetor com os custos garantidos avaliados no grid
%   .realCosts : Vetor com os custos reais avaliados no grid
%   .alpha     : Matriz com os pontos do simplex testados
%
% Date: 28/03/2026
% Author: glauco.gcps@gmail.com

%   Bloco de Validação de Entrada  
if nargin < 7
    error('h2_mjls_controle:NotEnoughInputs', 'São necessários A, B, E, C, D, Gamma e mu.');
end

input_cells = {A, B, E, C, D};
input_names = {'A', 'B', 'E', 'C', 'D'};
num_modos = numel(A); % sigma (Número de modos Markovianos)

% 1. Valida tipos e consistência do número de modos
for i = 1:length(input_cells)
    if ~iscell(input_cells{i})
        error('h2_mjls_controle:InvalidInputType', 'O argumento %s deve ser um cell array.', input_names{i});
    end
    if numel(input_cells{i}) ~= num_modos
        error('h2_mjls_controle:DimensionMismatch', 'Todos os cell arrays devem ter o mesmo número de modos (%d).', num_modos);
    end
end

if ~iscell(Gamma)
    error('h2_mjls_controle:InvalidGammaType', 'Gamma deve ser um cell array contendo os vértices do politopo.');
end
num_vertices = numel(Gamma); % N (Número de vértices do politopo)

% 2. Extrai dimensões baseando-se no primeiro modo
nx = size(A{1}, 1);
nu = size(B{1}, 2);
nw = size(E{1}, 2);
nz = size(C{1}, 1);

dims_A = size(A{1});
dims_B = size(B{1});
dims_E = size(E{1});
dims_C = size(C{1});
dims_D = size(D{1});

% 3. Validações de coerência das dimensões
if dims_A(1) ~= dims_A(2)
    error('h2_mjls_controle:InvalidADim', 'Matrizes A devem ser quadradas (nx x nx).');
end
if dims_B(1) ~= nx
    error('h2_mjls_controle:IncompatibleDims', 'O número de linhas de B deve ser igual ao de A (nx = %d).', nx);
end
if dims_E(1) ~= nx
    error('h2_mjls_controle:IncompatibleDims', 'O número de linhas de E deve ser igual ao de A (nx = %d).', nx);
end
if dims_C(2) ~= nx
    error('h2_mjls_controle:IncompatibleDims', 'O número de colunas de C deve ser igual ao de A (nx = %d).', nx);
end
if dims_D(1) ~= nz || dims_D(2) ~= nu
    error('h2_mjls_controle:IncompatibleDims', 'Dimensões de D devem ser (nz x nu) -> (%d x %d).', nz, nu);
end

% 4. Valida conteúdo numérico e constância das dimensões nos modos
for k = 1:num_modos
    if ~isnumeric(A{k}) || ~isnumeric(B{k}) || ~isnumeric(E{k}) || ~isnumeric(C{k}) || ~isnumeric(D{k})
        error('h2_mjls_controle:NonNumeric', 'O conteúdo de todas as células no modo %d deve ser numérico.', k);
    end
    if ~isequal(size(A{k}), dims_A) || ~isequal(size(B{k}), dims_B) || ...
       ~isequal(size(E{k}), dims_E) || ~isequal(size(C{k}), dims_C) || ~isequal(size(D{k}), dims_D)
        error('h2_mjls_controle:VaryingDims', 'As dimensões das matrizes não podem variar entre os modos (verifique o modo %d).', k);
    end
end

% 5. Validação das Propriedades Markovianas
for r = 1:num_vertices
    if size(Gamma{r}, 1) ~= num_modos || size(Gamma{r}, 2) ~= num_modos
        error('h2_mjls_controle:InvalidGamma', 'Cada vértice de Gamma deve ser (%d x %d).', num_modos, num_modos);
    end
    if any(Gamma{r}(:) < 0) || any(abs(sum(Gamma{r}, 2) - 1) > 1e-6)
        error('h2_mjls_controle:InvalidGammaProb', 'O vértice %d de Gamma não é uma matriz estocástica válida.', r);
    end
end

if length(mu) ~= num_modos || any(mu < 0) || abs(sum(mu) - 1) > 1e-6
    error('h2_mjls_controle:InvalidMu', 'mu deve ser uma distribuição de probabilidade válida.');
end

%   Tratamento das Opções (varargin)
options = [];
if nargin > 7
    if nargin == 8 && isstruct(varargin{1})
        options = varargin{1};
    else
        try
            options = struct(varargin{:});
        catch
            error('h2_mjls_controle:InvalidOptions', 'Opções inválidas. Use struct ou pares nome-valor.');
        end
    end
end

if ~isfield(options,'solver'), options.solver = 'sedumi'; end
if ~isfield(options,'deg'), options.deg = 1; end
if ~isfield(options,'h2'), options.h2 = 0; end
if ~isfield(options,'tolerance'), options.tolerance = 1e-7; end
if ~isfield(options,'op'), options.op = 1; end
if ~isfield(options,'degrho'), options.degrho = 0; end
if ~isfield(options,'indep'), options.indep = 0; end

%   Início da Lógica da Função
sigma = num_modos; 
N = num_vertices; 

% Construção da variável ROLMIP para a matriz de transição
Gamma_rol = cell(1, sigma);
for i=1:sigma
    Gaux = cell(1, N);
    for r=1:N
        % Transpor a coluna i para formar um bloco-linha
        Gaux{r} = kron(Gamma{r}(:,i)', eye(nx));
    end
    Gamma_rol{i} = rolmipvar(Gaux, strcat('Pr_alpha_',num2str(i)), N, 1);
end

output.cpusec_m = clock;
LMIs = [];

% Configuração da variável de desempenho (rho/gamma)
if options.h2 == 0
    if options.op == 1
        rho = sdpvar(1, 1);
        rho_alpha = rolmipvar(rho, 'rho', N, 0); 
        obj = rho;
    else        
        powers = gen_coefs(N, options.degrho);
        obj = [];
        rhos = cell(1, size(powers, 1));
        for i = 1:size(powers, 1)
            rhos{i} = sdpvar(1, 1);
            c_val = factorial(options.degrho) / prod(factorial(powers(i, :)));
            if isempty(obj)
                obj = rhos{i} * (1/c_val);
            else
                obj = rhos{i} * (1/c_val) + obj;
            end
        end    
        rho_alpha = rolmipvar(rhos, 'rho', N, options.degrho);
    end
else
    rho = options.h2 * options.h2;
    rho_alpha = rolmipvar(rho, 'rho', N, 0);
    options.tolerance = 0;
    obj = [];
end

% Variáveis de Lyapunov
Xcal = []; 
for i=1:sigma
    P{i} = rolmipvar(nx, nx, strcat('P',num2str(i)), 'symmetric', N, options.deg);
    if i == 1
        Xcal = P{i};
    else
        Xcal = [Xcal; P{i}];
    end
    LMIs = [LMIs, P{i} >= 0];
end

if (options.indep)
    H_in = rolmipvar(nx, nx, 'H_in', 'full', N, 0);
    Z_in = rolmipvar(nu, nx, 'Z_in', 'full', N, 0);
end

soma = 0;
for i=1:sigma
    W{i} = rolmipvar(nz, nz, strcat('W',num2str(i)), 'symmetric', N, options.deg);
    
    if (options.indep)
        H{i} = H_in;
        Z{i} = Z_in;
    else
        H{i} = rolmipvar(nx, nx, strcat('H',num2str(i)), 'full', N, 0);
        Z{i} = rolmipvar(nu, nx, strcat('Z',num2str(i)), 'full', N, 0);
    end
    
    soma = soma + trace(W{i});
    
    Q11 = W{i};
    Q12 = C{i}*H{i} + D{i}*Z{i};
    Q22 = -Gamma_rol{i}*Xcal + H{i} + H{i}';
    Q = [Q11, Q12;
         Q12', Q22];
     
    R11 = P{i} - mu(i)*E{i}*E{i}';
    R12 = A{i}*H{i} + B{i}*Z{i};
    R22 = -Gamma_rol{i}*Xcal + H{i} + H{i}';
    R = [R11, R12;
         R12', R22];
         
    LMIs = [LMIs, Q >= 0, R >= 0];
end 
LMIs = [LMIs, soma - rho_alpha <= 0];

output.cpusec_m = etime(clock, output.cpusec_m);
output.L = 0;
for i=1:size(LMIs,1)
    output.L = output.L + size(sdpvar(LMIs(i)),1);
end
output.V = size(getvariables(LMIs),2);

sol = solvesdp(LMIs, obj, sdpsettings('verbose', 0, 'solver', options.solver));
output.cpusec = sol.solvertime;
p = min(checkset(LMIs));
output.delta = p;

% Inicializa as saídas no padrão correto
output.feas = 0;
output.wc = Inf;
output.normSysWc = Inf;
output.alpha = [];
output.gcosts = [];
output.realCosts = [];
output.wc_alpha = [];
output.K = cell(1, sigma);

%   Validação (Grid Search em todo o Simplex)
if p > -options.tolerance
    output.feas = 1;
    
    % Extrai as matrizes P e calcula os ganhos do controlador K_i
    for i=1:sigma
        output.P{i} = double(P{i});
        H_val = double(H{i});
        Z_val = double(Z{i});
        output.K{i} = Z_val / H_val; % K_i = Z_i * H_i^-1
    end
    
    % Pré-calcula as matrizes do sistema em malha fechada (Fixas para todo alpha!)
    Acl = cell(1, sigma);
    Ecl = cell(1, sigma);
    Ccl = cell(1, sigma);
    Ezcl = cell(1, sigma);
    
    for i=1:sigma
        Acl{i} = A{i} + B{i}*output.K{i};
        Ecl{i} = E{i};
        Ccl{i} = C{i} + D{i}*output.K{i};
        Ezcl{i} = zeros(nz, nw); % Em H2 padrão, Dzw = 0
    end
    
    % Salva as matrizes de malha fechada na estrutura de saída
    output.Acl = Acl;
    output.Ecl = Ecl;
    output.Ccl = Ccl;
    output.Ezcl = Ezcl;
    
    % Variável ROLMIP para avaliar Gamma no grid
    Gamma_alpha = rolmipvar(Gamma, 'Ga', N, 1);
    
    passo_grid = 0.05; 
    pontos_alpha = particao_simplex(N, passo_grid);
    output.alpha = zeros(length(pontos_alpha), N);
    
    maxNorm = -1e10;
    wc_vec = [];
    
    for k = 1:length(pontos_alpha)
        alpha_vec = pontos_alpha{k}; 
        s = alpha_vec'; 
        output.alpha(k, :) = s;
        
        % 1. Custo Garantido (via LMI parametrizada)
        val = double(evalpar(rho_alpha, {s}));
        v_garantido = sqrt(max(0, val));
        output.gcosts = [output.gcosts; v_garantido];
        
        % Rastreia o pior caso garantido
        if v_garantido > maxNorm
            maxNorm = v_garantido;
            wc_vec = s;
        end
        
        % 2. Custo Real (Sistema em MALHA FECHADA Congelado)
        % Avalia a matriz de transição no ponto do simplex
        Gamma_aux = double(evalpar(Gamma_alpha, {s}));
        
        % Chama a função de análise nominal com as matrizes de malha fechada
%         out_nominal = h2_lmi_mjls_d(Acl, Ecl, Ccl, Ezcl, Gamma_aux, mu);      
%         if out_nominal.feas == 1
%             output.realCosts = [output.realCosts; out_nominal.h2];

        
        % Chama a função analítica via Kronecker (MUITO mais rápido)
        [h2_real, feas_real] = calc_h2_nominal_kronecker(Acl, Ecl, Ccl, Gamma_aux, mu);       
        if feas_real == 1
            output.realCosts = [output.realCosts; h2_real];

        else
            output.realCosts = [output.realCosts; NaN];
        end
    end
    
    output.wc = maxNorm; 
    output.wc_alpha = wc_vec;
    output.h2 = maxNorm;
    
    % Calcula a norma real do pior caso garantido
    Gamma_wc = double(evalpar(Gamma_alpha, {wc_vec}));
%     out_wc = h2_lmi_mjls_d(Acl, Ecl, Ccl, Ezcl, Gamma_wc, mu);
%     if out_wc.feas == 1
%         output.normSysWc = out_wc.h2;

    [h2_wc_real, feas_wc_real] = calc_h2_nominal_kronecker(Acl, Ecl, Ccl, Gamma_wc, mu);
    if feas_wc_real == 1
        output.normSysWc = h2_wc_real;
    else
        output.normSysWc = NaN;
    end
end
end

function [h2_norm, feas] = calc_h2_nominal_kronecker_primal(A, E, C, Gamma, mu)
% Calcula a norma H2 nominal de um MJLS discreto via Produto de Kronecker
% Resolve o sistema de equações de Lyapunov acopladas:
% P_i - A_i^T * sum(Gamma_ij * P_j) * A_i = C_i^T * C_i

    sigma = length(A);
    nx = size(A{1}, 1);
    nx2 = nx^2;
    
    % Montagem do grande sistema linear M * X = V
    M = zeros(sigma * nx2, sigma * nx2);
    V = zeros(sigma * nx2, 1);
    
    for i = 1:sigma
        % Termo independente: vec(C_i^T * C_i)
        V((i-1)*nx2 + 1 : i*nx2) = C{i}' * C{i}; % O MATLAB vetoriza automaticamente matrizes em atribuições de coluna
        
        % Bloco diagonal principal (Identidade)
        M((i-1)*nx2 + 1 : i*nx2, (i-1)*nx2 + 1 : i*nx2) = eye(nx2);
        
        % Termo acoplado: (A_i^T \otimes A_i^T) * Gamma_ij
        Ai_kron = kron(A{i}', A{i}');
        for j = 1:sigma
            M((i-1)*nx2 + 1 : i*nx2, (j-1)*nx2 + 1 : j*nx2) = ...
                M((i-1)*nx2 + 1 : i*nx2, (j-1)*nx2 + 1 : j*nx2) - Gamma(i,j) * Ai_kron;
        end
    end
    
    % Verifica estabilidade (se a matriz M é singular ou mal condicionada)
    if rcond(M) < 1e-15
        h2_norm = NaN;
        feas = 0;
        return;
    end
    
    % Resolve o sistema linear X = M \ V
    X = M \ V;
    
    % Reconstrói as matrizes P_i e calcula a norma H2
    h2_sq = 0;
    for i = 1:sigma
        % Extrai o vec(P_i) e remodela para matriz nx x nx
        Pi_vec = X((i-1)*nx2 + 1 : i*nx2);
        Pi = reshape(Pi_vec, nx, nx);
        
        % Garante simetria (erros numéricos)
        Pi = (Pi + Pi') / 2;
        
        % Soma o custo: tr(E_i^T * P_i * E_i * mu_i)
        h2_sq = h2_sq + trace(E{i}' * Pi * E{i}) * mu(i);
    end
    
    % Se o custo for negativo (erro numérico em sistemas instáveis), é inviável
    if h2_sq < 0
        h2_norm = NaN;
        feas = 0;
    else
        h2_norm = sqrt(h2_sq);
        feas = 1;
    end
end

function [h2_norm, feas] = calc_h2_nominal_kronecker(A, E, C, Gamma, mu)
% Calcula a norma H2 nominal de um MJLS discreto via Produto de Kronecker
% Resolve o sistema de equações de Lyapunov acopladas (Formulação Dual/Covariância):
% P_i - A_i * sum(Gamma_ji * P_j) * A_i^T = mu_i * E_i * E_i^T

    sigma = length(A);
    nx = size(A{1}, 1);
    nx2 = nx^2;
    
    % Montagem do grande sistema linear M * X = V
    M = zeros(sigma * nx2, sigma * nx2);
    V = zeros(sigma * nx2, 1);
    
    for i = 1:sigma
        % Termo independente: vec(mu_i * E_i * E_i^T)
        Ei_term = mu(i) * (E{i} * E{i}');
        V((i-1)*nx2 + 1 : i*nx2) = Ei_term(:); % O (:) vetoriza a matriz
        
        % Bloco diagonal principal (Identidade)
        M((i-1)*nx2 + 1 : i*nx2, (i-1)*nx2 + 1 : i*nx2) = eye(nx2);
        
        % Termo acoplado: (A_i \otimes A_i) * Gamma_ji
        % usamos A_i sem transpor, pois é a formulação de covariância
        Ai_kron = kron(A{i}, A{i});
        for j = 1:sigma
            % a soma é em j, e o coeficiente é Gamma(j,i)
            M((i-1)*nx2 + 1 : i*nx2, (j-1)*nx2 + 1 : j*nx2) = ...
                M((i-1)*nx2 + 1 : i*nx2, (j-1)*nx2 + 1 : j*nx2) - Gamma(j,i) * Ai_kron;
        end
    end
    
    % Verifica estabilidade (se a matriz M é singular ou mal condicionada)
    if rcond(M) < 1e-15
        h2_norm = NaN;
        feas = 0;
        return;
    end
    
    % Resolve o sistema linear X = M \ V
    X = M \ V;
    
    % Reconstrói as matrizes P_i
    P = cell(1, sigma);
    for i = 1:sigma
        Pi_vec = X((i-1)*nx2 + 1 : i*nx2);
        Pi = reshape(Pi_vec, nx, nx);
        P{i} = (Pi + Pi') / 2; % Garante simetria contra erros numéricos
    end
    
    % Calcula a covariância real S_i = sum_{j=1}^sigma Gamma(j,i) * P_j
    S = cell(1, sigma);
    for i = 1:sigma
        S{i} = zeros(nx, nx);
        for j = 1:sigma
            S{i} = S{i} + Gamma(j,i) * P{j};
        end
    end
    
    % Calcula o custo H2 final: sum_{i=1}^sigma trace(C_i * S_i * C_i^T)
    h2_sq = 0;
    for i = 1:sigma
        h2_sq = h2_sq + trace(C{i} * S{i} * C{i}');
    end
    
    % Se o custo for negativo (erro numérico em sistemas instáveis), é inviável
    if h2_sq < 0
        h2_norm = NaN;
        feas = 0;
    else
        h2_norm = sqrt(h2_sq);
        feas = 1;
    end
end
