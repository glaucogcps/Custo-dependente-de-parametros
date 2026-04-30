function output = hinf_lmi_mjls_d_incerto_analise(A,E,Cz,Ez,Gamma,mu,varargin)
% function output = hinf_lmi_mjls_d_incerto_analise(A,E,Cz,Ez,Gamma,mu,varargin)
%
% Evaluate the H-infinity norm of a discrete-time MJLS with polytopic uncertainties. 
% Based on Lemma 1 for uncertain transition probabilities.
% The LMIs are programmed using YALMIP/ROLMIP and can be solved by any 
% LMI solver supported by YALMIP (SeDuMi is the default). 
%
% input:  (A,E,Cz,Ez) -> state-space matrices (cell arrays for modes)
%         Gamma       -> cell array with vertices of the probability matrix polytope
%         mu          -> initial probability distribution
%
% varargin (Opções):
%   .solver    : Solver LMI a ser utilizado (default: 'sedumi')
%   .deg       : Grau do polinômio da matriz de Lyapunov P (default: 1)
%   .hinf      : Limite superior desejado para a norma H-inf. Se for 0, a 
%                função minimiza a norma H-inf (default: 0)
%   .tolerance : Tolerância numérica para a verificação de viabilidade (default: 1e-7)
%   .op        : Tipo de otimização (0 = integral/média, 1 = pico) (default: 1)
%   .degrho    : Grau do polinômio da variável de desempenho rho (default: 0)
%        
% output: 
%   .feas      : Status de viabilidade (1 = viável, 0 = inviável)
%   .hinf      : H-inf norm guaranteed worst-case (0 if unfeasible)
%   .wc        : Pior custo H-inf garantido encontrado no grid
%   .normSysWc : Custo H-inf real do sistema no pior caso
%   .cpusec    : cpu time to solve the LMIs (seconds)
%   .cpusec_m  : cpu time to mount the LMIs (seconds)
%   .P         : Lyapunov matrix (ROLMIP variables)
%   .V         : number of decision variables
%   .L         : number of LMIs rows
%   .gcosts    : Vetor com os custos garantidos avaliados no grid
%   .realCosts : Vetor com os custos reais avaliados no grid
%   .alpha     : Matriz com os pontos do simplex testados
%
% Date: 04/04/2026
% Author: glauco.gcps@gmail.com

if nargin < 6
    error('hinf_mjls_incerto:NotEnoughInputs', 'São necessários A, E, Cz, Ez, Gamma e mu.');
end

input_cells = {A, E, Cz, Ez};
input_names = {'A', 'E', 'Cz', 'Ez'};
num_modos = numel(A); % sigma (Número de modos Markovianos)

% 1. Valida tipos e consistência do número de modos para matrizes do sistema
for i = 1:length(input_cells)
    if ~iscell(input_cells{i})
        error('hinf_mjls_incerto:InvalidInputType', 'O argumento %s deve ser um cell array.', input_names{i});
    end
    if numel(input_cells{i}) ~= num_modos
        error('hinf_mjls_incerto:DimensionMismatch', 'Todos os cell arrays do sistema devem ter o mesmo número de modos (%d).', num_modos);
    end
end

if ~iscell(Gamma)
    error('hinf_mjls_incerto:InvalidGammaType', 'Gamma deve ser um cell array contendo os vértices do politopo de transição.');
end

num_vertices = numel(Gamma); % N (Número de vértices do politopo)

% 2. Extrai dimensões baseando-se no primeiro modo
nx = size(A{1}, 1);
nw = size(E{1}, 2);
nz = size(Cz{1}, 1);

dims_A  = size(A{1});
dims_E  = size(E{1});
dims_Cz = size(Cz{1});
dims_Ez = size(Ez{1});

% 3. Validações de coerência das dimensões
if dims_A(1) ~= dims_A(2)
    error('hinf_mjls_incerto:InvalidADim', 'Matrizes A devem ser quadradas (nx x nx).');
end
if dims_E(1) ~= nx
    error('hinf_mjls_incerto:IncompatibleDims', 'O número de linhas de E deve ser igual ao de A (nx = %d).', nx);
end
if dims_Cz(2) ~= nx
    error('hinf_mjls_incerto:IncompatibleDims', 'O número de colunas de Cz deve ser igual ao de A (nx = %d).', nx);
end
if dims_Ez(1) ~= nz || dims_Ez(2) ~= nw
    error('hinf_mjls_incerto:IncompatibleDims', 'Dimensões de Ez devem ser (nz x nw) -> (%d x %d).', nz, nw);
end

% 4. Valida conteúdo numérico e constância das dimensões nos modos
for k = 1:num_modos
    if ~isnumeric(A{k}) || ~isnumeric(E{k}) || ~isnumeric(Cz{k}) || ~isnumeric(Ez{k})
        error('hinf_mjls_incerto:NonNumeric', 'O conteúdo de todas as células no modo %d deve ser numérico.', k);
    end
    if ~isequal(size(A{k}), dims_A) || ~isequal(size(E{k}), dims_E) || ...
       ~isequal(size(Cz{k}), dims_Cz) || ~isequal(size(Ez{k}), dims_Ez)
        error('hinf_mjls_incerto:VaryingDims', 'As dimensões das matrizes não podem variar entre os modos (verifique o modo %d).', k);
    end
end

% 5. Validação das Propriedades Markovianas nos Vértices (Gamma) e mu
for r = 1:num_vertices
    if size(Gamma{r}, 1) ~= num_modos || size(Gamma{r}, 2) ~= num_modos
        error('hinf_mjls_incerto:InvalidGamma', 'Cada vértice de Gamma deve ser quadrado com dimensão igual ao número de modos (%d x %d).', num_modos, num_modos);
    end
    if any(Gamma{r}(:) < 0) || any(abs(sum(Gamma{r}, 2) - 1) > 1e-6)
        error('hinf_mjls_incerto:InvalidGammaProb', 'O vértice %d de Gamma não é uma matriz estocástica válida.', r);
    end
end

if length(mu) ~= num_modos
    error('hinf_mjls_incerto:InvalidMu', 'O vetor mu deve ter comprimento igual ao número de modos (%d).', num_modos);
end
if any(mu < 0) || abs(sum(mu) - 1) > 1e-6
    error('hinf_mjls_incerto:InvalidMuProb', 'mu deve ser uma distribuição de probabilidade válida (elementos >= 0 e soma = 1).');
end

%   Tratamento das Opções (varargin)
options = [];
if nargin > 6
    if nargin == 7 && isstruct(varargin{1})
        options = varargin{1};
    else
        try
            options = struct(varargin{:});
        catch
            error('hinf_mjls_incerto:InvalidOptions', 'Opções inválidas. Use struct ou pares nome-valor.');
        end
    end
end

if ~isfield(options,'solver'), options.solver = 'sedumi'; end
if ~isfield(options,'deg'), options.deg = 1; end
if ~isfield(options,'hinf'), options.hinf = 0; end % Corrigido de h2 para hinf
if ~isfield(options,'tolerance'), options.tolerance = 1e-7; end
if ~isfield(options,'op'), options.op = 1; end
if ~isfield(options,'degrho'), options.degrho = 0; end

%   Início da Lógica da Função
sigma = size(Gamma{1},1); %operation modes
N = size(Gamma,2); %vertices

output.cpusec_m = clock;

%new LMI system
LMIs = [];

% Configuração da variável de desempenho (rho/gamma)
if options.hinf > 0
    % Apenas checagem de viabilidade para um H-inf dado
    rho_alpha = options.hinf * options.hinf;
    options.tolerance = 0;
    obj = [];
else
    % Minimização da norma
    if options.op == 1
        % Caso op=1: rho é uma variável escalar constante (pico)
        rho = sdpvar(1, 1);
        rho_alpha = rolmipvar(rho, 'rho', N, 0); 
        obj = rho;
    else        
        % Caso op=0: rho é um polinômio rho(alpha)
        powers = gen_coefs(N, options.degrho);
        obj = [];
        rhos = cell(1, size(powers, 1));
        
        for i = 1:size(powers, 1)
            rhos{i} = sdpvar(1, 1);
            % Coeficiente binomial para ponderação na função objetivo (integral/média)
            c_val = factorial(options.degrho) / prod(factorial(powers(i, :)));
            
            if isempty(obj)
                obj = rhos{i} * (1/c_val);
            else
                obj = rhos{i} * (1/c_val) + obj;
            end
        end    
        rho_alpha = rolmipvar(rhos, 'rho', N, options.degrho);
    end
end

Xcal = []; 
for i=1:sigma
    P{i} = rolmipvar(nx,nx,strcat('P',num2str(i)),'symmetric', N, options.deg);
    Xcal = [Xcal; P{i}];
    LMIs = [LMIs, P{i} >= 0];
end

% Cálculo da esperança matemática vetorizado com Kronecker
for i=1:sigma
    for r=1:N
        Gaux{r} = kron(Gamma{r}(i,:),eye(nx));
    end
    Lambda{i} = rolmipvar(Gaux,strcat('Pr_alpha_',num2str(i)),N,1);
    Pp{i} = Lambda{i} * Xcal;
end

for i=1:sigma
    Z11 = A{i};
    Z12 = E{i};
    Z21 = Cz{i};
    Z22 = Ez{i};
    Z = [Z11, Z12;
         Z21, Z22];
    
    X = Z';
    
    Y11 = Pp{i};
    Y12 = zeros(nx, nz);
    Y22 = eye(nz);
    Y = [Y11, Y12;
         Y12', Y22];
    
    W11 = P{i};
    W12 = zeros(nx, nw);
    W21 = zeros(nw, nx);
    W22 = rho_alpha * eye(nw);
    W = [W11, W12;
         W21, W22];
         
    LMIs = [LMIs, X*Y*Z - W <= 0];
end 

output.cpusec_m = etime(clock,output.cpusec_m);

% determine the number of LMI rows
output.L = 0;
for i=1:size(LMIs,1)
    output.L = output.L + size(sdpvar(LMIs(i)),1);
end

% determine the number of scalar variables
output.V = size(getvariables(LMIs),2);

sol = solvesdp(LMIs, obj, sdpsettings('verbose', 0, 'solver', options.solver));
output.cpusec = sol.solvertime;

p = min(checkset(LMIs));
output.delta = p;

% Inicializa as saídas no padrão correto e separado
output.feas = 0;
output.hinf = 0;
output.wc = Inf;
output.normSysWc = Inf;
output.alpha = [];
output.gcosts = [];
output.realCosts = [];
output.wc_alpha = [];

%   Validação (Grid Search em todo o Simplex)
if p > -options.tolerance
    output.feas = 1;
    
    % Extrai as matrizes P (células do ROLMIP)
    for i=1:sigma
        output.P{i} = double(P{i});
    end
    
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
        
        % 1. Custo Garantido pela LMI parametrizada
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
        
        % 2. Custo Real usando a função de análise nominal
        Gamma_aux = double(evalpar(Gamma_alpha, {s}));
        out_nominal = hinf_lmi_mjls_d(A, E, Cz, Ez, Gamma_aux, mu);
        
        % Preenche o vetor realCosts (separado)
        if out_nominal.feas == 1
            output.realCosts = [output.realCosts; out_nominal.hinf]; 
        else
            output.realCosts = [output.realCosts; NaN];
        end
    end
    
    output.wc = maxNorm; % Pior caso garantido no grid
    output.wc_alpha = wc_vec;
    output.hinf = maxNorm; 
    
    % Calcula a norma real do pior caso garantido
    Gamma_wc = double(evalpar(Gamma_alpha, {wc_vec}));
    out_wc = hinf_lmi_mjls_d(A, E, Cz, Ez, Gamma_wc, mu);
    if out_wc.feas == 1
        output.normSysWc = out_wc.hinf; 
    else
        output.normSysWc = NaN;
    end
end
end