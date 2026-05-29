function output = hinf_lmi_mjls_d_incerto_controle(A, B, E, Cz, Dz, Ez, Gamma, mu, varargin)
% function output = hinf_lmi_mjls_d_incerto_controle(A, B, E, Cz, Dz, Ez, Gamma, mu, varargin)
%
% Synthesize an H-infinity guaranteed cost controller (u = K_i x) for a discrete-time 
% MJLS with uncertain probability matrix (polytopic).
%
% input:  (A,B,E,Cz,Dz,Ez) -> state-space matrices (cell arrays for modes)
%         Gamma            -> cell array with vertices of the probability matrix polytope
%         mu               -> initial probability distribution
%
% varargin (Opções):
%   .solver    : Solver LMI a ser utilizado (default: 'sedumi')
%   .deg       : Grau do polinômio da matriz de Lyapunov P (default: 1)
%   .hinf      : Limite superior desejado para a norma H-inf. Se for 0, a 
%                função minimiza a norma H-inf (default: 0)
%   .tolerance : Tolerância numérica para a verificação de viabilidade (default: 1e-7)
%   .op        : Tipo de otimização (0 = integral/média, 1 = pico) (default: 1)
%   .degrho    : Grau do polinômio da variável de desempenho rho (default: 0)
%   .indep     : Controlador independente de modo (1 = sim, 0 = não) (default: 0)
%        
% Saída (output):
%   .feas      : Status de viabilidade (1 = viável, 0 = inviável)
%   .hinf      : Pior custo H-inf garantido (mantido para retrocompatibilidade)
%   .wc        : Pior custo H-inf garantido encontrado no grid
%   .normSysWc : Custo H-inf real do sistema em malha fechada no pior caso
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
% Date: 04/04/2026
% Author: glauco.gcps@gmail.com

%   Bloco de Validação de Entrada  
if nargin < 8
    error('hinf_mjls_controle:NotEnoughInputs', 'São necessários A, B, E, Cz, Dz, Ez, Gamma e mu.');
end

input_cells = {A, B, E, Cz, Dz, Ez};
input_names = {'A', 'B', 'E', 'Cz', 'Dz', 'Ez'};
num_modos = numel(A); % sigma (Número de modos Markovianos)

% 1. Valida tipos e consistência do número de modos
for i = 1:length(input_cells)
    if ~iscell(input_cells{i})
        error('hinf_mjls_controle:InvalidInputType', 'O argumento %s deve ser um cell array.', input_names{i});
    end
    if numel(input_cells{i}) ~= num_modos
        error('hinf_mjls_controle:DimensionMismatch', 'Todos os cell arrays devem ter o mesmo número de modos (%d).', num_modos);
    end
end

if ~iscell(Gamma)
    error('hinf_mjls_controle:InvalidGammaType', 'Gamma deve ser um cell array contendo os vértices do politopo.');
end
num_vertices = numel(Gamma); % N (Número de vértices do politopo)

% 2. Extrai dimensões baseando-se no primeiro modo
nx = size(A{1}, 1);
nu = size(B{1}, 2);
nw = size(E{1}, 2);
nz = size(Cz{1}, 1);

dims_A  = size(A{1});
dims_B  = size(B{1});
dims_E  = size(E{1});
dims_Cz = size(Cz{1});
dims_Dz = size(Dz{1});
dims_Ez = size(Ez{1});

% 3. Validações de coerência das dimensões
if dims_A(1) ~= dims_A(2)
    error('hinf_mjls_controle:InvalidADim', 'Matrizes A devem ser quadradas (nx x nx).');
end
if dims_B(1) ~= nx
    error('hinf_mjls_controle:IncompatibleDims', 'O número de linhas de B deve ser igual ao de A (nx = %d).', nx);
end
if dims_E(1) ~= nx
    error('hinf_mjls_controle:IncompatibleDims', 'O número de linhas de E deve ser igual ao de A (nx = %d).', nx);
end
if dims_Cz(2) ~= nx
    error('hinf_mjls_controle:IncompatibleDims', 'O número de colunas de Cz deve ser igual ao de A (nx = %d).', nx);
end
if dims_Dz(1) ~= nz || dims_Dz(2) ~= nu
    error('hinf_mjls_controle:IncompatibleDims', 'Dimensões de Dz devem ser (nz x nu) -> (%d x %d).', nz, nu);
end
if dims_Ez(1) ~= nz || dims_Ez(2) ~= nw
    error('hinf_mjls_controle:IncompatibleDims', 'Dimensões de Ez devem ser (nz x nw) -> (%d x %d).', nz, nw);
end

% 4. Valida conteúdo numérico e constância das dimensões nos modos
for k = 1:num_modos
    if ~isnumeric(A{k}) || ~isnumeric(B{k}) || ~isnumeric(E{k}) || ~isnumeric(Cz{k}) || ~isnumeric(Dz{k}) || ~isnumeric(Ez{k})
        error('hinf_mjls_controle:NonNumeric', 'O conteúdo de todas as células no modo %d deve ser numérico.', k);
    end
    if ~isequal(size(A{k}), dims_A) || ~isequal(size(B{k}), dims_B) || ...
       ~isequal(size(E{k}), dims_E) || ~isequal(size(Cz{k}), dims_Cz) || ...
       ~isequal(size(Dz{k}), dims_Dz) || ~isequal(size(Ez{k}), dims_Ez)
        error('hinf_mjls_controle:VaryingDims', 'As dimensões das matrizes não podem variar entre os modos (verifique o modo %d).', k);
    end
end

% 5. Validação das Propriedades Markovianas
for r = 1:num_vertices
    if size(Gamma{r}, 1) ~= num_modos || size(Gamma{r}, 2) ~= num_modos
        error('hinf_mjls_controle:InvalidGamma', 'Cada vértice de Gamma deve ser (%d x %d).', num_modos, num_modos);
    end
    if any(Gamma{r}(:) < 0) || any(abs(sum(Gamma{r}, 2) - 1) > 1e-6)
        error('hinf_mjls_controle:InvalidGammaProb', 'O vértice %d de Gamma não é uma matriz estocástica válida.', r);
    end
end

if length(mu) ~= num_modos || any(mu < 0) || abs(sum(mu) - 1) > 1e-6
    error('hinf_mjls_controle:InvalidMu', 'mu deve ser uma distribuição de probabilidade válida.');
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
            error('hinf_mjls_controle:InvalidOptions', 'Opções inválidas. Use struct ou pares nome-valor.');
        end
    end
end

if ~isfield(options,'solver'), options.solver = 'sedumi'; end
if ~isfield(options,'deg'), options.deg = 1; end
if ~isfield(options,'hinf'), options.hinf = 0; end
if ~isfield(options,'tolerance'), options.tolerance = 1e-7; end
if ~isfield(options,'op'), options.op = 1; end
if ~isfield(options,'degrho'), options.degrho = 0; end
if ~isfield(options,'indep'), options.indep = 0; end

%   Início da Lógica da Função
sigma = num_modos; 
N = num_vertices; 

% Construção da variável ROLMIP para a matriz de transição
% Gamma_rol = cell(1, sigma);
% Gamma_diag = [];
% for i=1:sigma
%     Gaux = cell(1, N);
%     for r=1:N
%         Gaux{r} = kron(Gamma{r}(i,:), eye(nx));
%     end
%     Gamma_rol{i} = rolmipvar(Gaux, strcat('Pr_alpha_',num2str(i)), N, 1);
%     Gamma_diag = blkdiag(Gamma_diag, Gamma_rol{i});
% end

Gamma_rol = cell(1, sigma);
Gamma_diag = cell(1, sigma);
for i=1:sigma
    Gaux = cell(1, N);
    Gdiag = cell(1, N);
    for r=1:N
        Gaux{r} = kron(Gamma{r}(i,:), eye(nx));
        Gauxd = [];
        for j = 1:sigma
            Gauxd = blkdiag(Gauxd, Gamma{r}(i,j)*eye(nx));
        end
        Gdiag{r} = Gauxd;
    end
    Gamma_rol{i} = rolmipvar(Gaux, strcat('Pr_alpha_',num2str(i)), N, 1);
    Gamma_diag{i} = rolmipvar(Gdiag, strcat('Pr_alpha_d',num2str(i)), N, 1);
end

output.cpusec_m = clock;
LMIs = [];

% Configuração da variável de desempenho (rho/gamma)
if options.hinf == 0
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
    rho = options.hinf * options.hinf;
    rho_alpha = rolmipvar(rho, 'rho', N, 0);
    options.tolerance = 0;
    obj = [];
end

if (options.indep)
    H_in = rolmipvar(nx, nx, 'H_in', 'full', N, 0);
    Z_in = rolmipvar(nu, nx, 'Z_in', 'full', N, 0);
end

% Variáveis de Lyapunov
Pbar = [];
for i=1:sigma
    P{i} = rolmipvar(nx, nx, strcat('P',num2str(i)), 'symmetric', N, options.deg);
    
    if (options.indep)
        H{i} = H_in;
        Z{i} = Z_in;
    else
        H{i} = rolmipvar(nx, nx, strcat('H',num2str(i)), 'full', N, 0);
        Z{i} = rolmipvar(nu, nx, strcat('Z',num2str(i)), 'full', N, 0);
    end
    
    if i == 1
        Pbar = P{i};
    else
        Pbar = blkdiag(Pbar, P{i});
    end
    LMIs = [LMIs, P{i} >= 0];
end

for i=1:sigma
    M11 = P{i} - H{i} - H{i}';
    M21 = zeros(nw, nx);
    M22 = -rho_alpha*eye(nw);
    M31 = Gamma_rol{i}'*(A{i}*H{i} + B{i}*Z{i});
    M32 = Gamma_rol{i}'*E{i};
    M33 = -Gamma_diag{i}*Pbar;
    M41 = Cz{i}*H{i} + Dz{i}*Z{i};
    M42 = Ez{i};
    M43 = zeros(nz, sigma*nx);
    M44 = -eye(nz);
    M = [M11, M21', M31', M41';
         M21, M22, M32', M42';
         M31, M32, M33, M43';
         M41, M42, M43, M44];
         
    LMIs = [LMIs, M <= 0];
end 


output.cpusec_m = etime(clock, output.cpusec_m);

output.L = 0;
for i=1:size(LMIs,1)
    output.L = output.L + size(sdpvar(LMIs(i)),1);
end
output.V = size(getvariables(LMIs),2);

sol = solvesdp(LMIs, obj, sdpsettings('verbose', 0, 'solver', options.solver));
output.cpusec = sol.solvertime;

p = min(checkset(LMIs));
output.delta = 0;
output.feas = 0;
output.hinf = 0;
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
        if sum(isnan(output.K{i}))
            output.K{i} = zeros(nu, nx);
        end
        Acl{i} = A{i} + B{i}*output.K{i};
        Ecl{i} = E{i};
        Ccl{i} = Cz{i} + Dz{i}*output.K{i};
        Ezcl{i} = Ez{i}; % Em H-infinito, Ez pode ser não-nulo
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
    
    % Pré-aloca a matriz alpha
    output.alpha = zeros(length(pontos_alpha), N);
    
    maxNorm = -1e10;
    wc_vec = [];
    
    for k = 1:length(pontos_alpha)
        alpha_vec = pontos_alpha{k}; % Vetor coluna
        s = alpha_vec'; % Vetor linha para o ROLMIP
        
        % Preenche a matriz alpha (separada)
        output.alpha(k, :) = s;
        
        % 1. Custo Garantido (via LMI parametrizada)
        if options.hinf > 0
            val = options.hinf^2;
        else
            val = double(evalpar(rho_alpha, {s}));
        end
        v_garantido = double(sqrt(max(0, val)));
        
        % Preenche o vetor gcosts (separado)
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
        out_nominal = hinf_lmi_mjls_d(Acl, Ecl, Ccl, Ezcl, Gamma_aux, mu);
        
        % Preenche o vetor realCosts (separado)
        if out_nominal.feas == 1
            output.realCosts = [output.realCosts; out_nominal.hinf];
        else
            output.realCosts = [output.realCosts; NaN];
        end
    end
    
    output.wc = maxNorm; 
    output.wc_alpha = wc_vec;
    output.hinf = maxNorm;
    
    % Calcula a norma real do pior caso garantido
    Gamma_wc = double(evalpar(Gamma_alpha, {wc_vec}));
    out_wc = hinf_lmi_mjls_d(Acl, Ecl, Ccl, Ezcl, Gamma_wc, mu);
    if out_wc.feas == 1
        output.normSysWc = out_wc.hinf;
    else
        output.normSysWc = NaN;
    end
end
end





% % Variáveis de Lyapunov
% Xcal = []; 
% for i=1:sigma
%     P{i} = rolmipvar(nx, nx, strcat('P',num2str(i)), 'symmetric', N, options.deg);
%     Xcal = [Xcal; P{i}]; % Empilha P_i verticalmente para formar Xcal
%     LMIs = [LMIs, P{i} >= 0];
% end
% 
% % Construção do operador (Lambda_i)
% for i=1:sigma
%     for r=1:N
%         Gaux{r} = kron(Gamma{r}(i,:), eye(nx));
%     end
%     Lambda{i} = rolmipvar(Gaux, strcat('Pr_alpha_',num2str(i)), N, 1);
% end
% 
% for i=1:sigma
%     if (options.indep)
%         H{i} = H_in;
%         Z{i} = Z_in;
%     else
%         H{i} = rolmipvar(nx, nx, strcat('H',num2str(i)), 'full', N, 0);
%         Z{i} = rolmipvar(nu, nx, strcat('Z',num2str(i)), 'full', N, 0);
%     end
%     
%     % Montagem exata do Teorema 2 (Forma Transposta para bater com seu código)
%     M11 = P{i} - H{i} - H{i}';
%     M12 = zeros(nx, nw);
%     
%     % O truque do Teorema 2: empilhar A_cl e E para todos os modos futuros
%     % Em vez de multiplicar por Gamma_rol, nós multiplicamos a matriz Xcal por Lambda
%     % Isso é matematicamente equivalente a usar Upsilon_i e X_i do artigo
%     
%     % Matriz Z do Teorema 2 (Dinâmica)
%     Z11 = A{i}*H{i} + B{i}*Z{i};
%     Z12 = E{i};
%     Z21 = Cz{i}*H{i} + Dz{i}*Z{i};
%     Z22 = Ez{i};
%     
%     Z_mat = [Z11, Z12;
%              Z21, Z22];
%              
%     % Matriz Y do Teorema 2 
%     Y11 = Lambda{i} * Xcal; % Isso faz sum(p_ij * P_j)
%     Y12 = zeros(nx, nz);
%     Y22 = eye(nz);
%     
%     Y_mat = [Y11, Y12;
%              Y12', Y22];
%              
%     % Matriz W do Teorema 2 (Lyapunov atual e Desempenho)
%     W11 = P{i} - H{i} - H{i}';
%     W12 = zeros(nx, nw);
%     W22 = -rho_alpha * eye(nw);
%     
%     W_mat = [W11, W12;
%              W12', W22];
%              
%     % A LMI final: Z^T * Y * Z + W < 0
%     % Pelo Complemento de Schur, isso é equivalente à LMI do Teorema 2
%     % E é muito mais eficiente numericamente no YALMIP
%     
%     % Montagem via Complemento de Schur (Forma expandida do Teorema 2)
%     LMI_i = [W11, W12, Z11', Z21';
%              W12', W22, Z12', Z22';
%              Z11, Z12, -Y11, Y12;
%              Z21, Z22, Y12', -Y22];
%              
%     LMIs = [LMIs, LMI_i <= 0];
% end