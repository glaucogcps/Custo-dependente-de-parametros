function out = hinf_lmi_d_incerto_controle(A, B, E, Cz, Dz, Ez, varargin)
% function out = hinf_lmi_d_incerto_controle(A, B, E, Cz, Dz, Ez, varargin)
%
% Sintetiza um controlador H infinito garantido (u = Kx) para um SISTEMA 
% LINEAR DISCRETO no tempo com incerteza politópica.
%
% Esta função implementa uma formulação relaxada baseada em multiplicadores
% onde a matriz de Lyapunov (W) e o custo (rho_alpha) são dependentes de 
% parâmetros (grau 'deg'), enquanto as variáveis de síntese (X e Z) são 
% de grau 0 (controlador robusto fixo para todo o politopo).
%
% O sistema em malha aberta é dado por:
% x(k+1) = A(alpha)x(k) + B(alpha)u(k) + E(alpha)w(k)
% z(k)   = Cz(alpha)x(k) + Dz(alpha)u(k) + Ez(alpha)w(k)
%
% Entradas (Cell arrays com N matrizes cada):
%   A  - Matriz de estado (nx x nx)
%   B  - Matriz de entrada de controle (nx x nu)
%   E  - Matriz de entrada de distúrbio (nx x nw)
%   Cz - Matriz de saída de desempenho (nz x nx)
%   Dz - Matriz de transmissão direta do controle para saída (nz x nu)
%   Ez - Matriz de transmissão direta do distúrbio para saída (nz x nw)
%
% varargin (Opções):
%   .solver   : Solver a ser utilizado (default: 'mosek')
%   .deg      : Grau do polinômio da matriz W (default: 1)
%   .degGamma : Grau do polinômio da variável de desempenho (default: 0)
%   .op       : Tipo de otimização (0 = integral/média, 1 = pico) (default: 0)
%   .verbose  : Nível de detalhe do solver (default: 0)
%
% Saída (out):
%   .feas      : Status de viabilidade (1 = viável, 0 = inviável)
%   .wc        : Pior custo Hinf garantido encontrado no grid
%   .normSysWc : Custo Hinf real do sistema em malha fechada no pior caso
%   .K_wc      : Ganho do controlador K = Z/X avaliado no pior caso
%   .Acl_wc... : Matrizes do sistema em malha fechada no pior caso
%   .gcosts    : Vetor com os custos garantidos avaliados no grid
%   .realCosts : Vetor com os custos reais avaliados no grid
%   .alpha     : Matriz com os pontos do simplex testados
%   .Z_cells, .X_cells : Variáveis de decisão para reconstruir K
%
% Date: 22/03/2026

%   Bloco de Validação de Entrada  
if nargin < 6
    error('hinf_controle:NotEnoughInputs', 'São necessários A, B, E, Cz, Dz, Ez (cell arrays de matrizes).');
end

input_cells = {A, B, E, Cz, Dz, Ez};
input_names = {'A', 'B', 'E', 'Cz', 'Dz', 'Ez'};
num_vertices = numel(A); % N (Número de vértices do politopo)

% 1. Valida tipos e consistência do número de vértices
for i = 1:length(input_cells)
    if ~iscell(input_cells{i})
        error('hinf_controle:InvalidInputType', 'O argumento %s deve ser um cell array.', input_names{i});
    end
    if numel(input_cells{i}) ~= num_vertices
        error('hinf_controle:DimensionMismatch', 'Todos os cell arrays devem ter o mesmo número de vértices (%d).', num_vertices);
    end
end

% 2. Extrai dimensões baseando-se no primeiro vértice
nx = size(A{1}, 1); % Número de estados
nu = size(B{1}, 2); % Número de entradas de controle
nw = size(E{1}, 2); % Número de entradas de distúrbio
nz = size(Cz{1}, 1); % Número de saídas de desempenho

dims_A  = size(A{1});
dims_B  = size(B{1});
dims_E  = size(E{1});
dims_Cz = size(Cz{1});
dims_Dz = size(Dz{1});
dims_Ez = size(Ez{1});

% 3. Validações de coerência das dimensões (Álgebra Linear)
if dims_A(1) ~= dims_A(2)
    error('hinf_controle:InvalidADim', 'Matrizes A devem ser quadradas (nx x nx).');
end
if dims_B(1) ~= nx
    error('hinf_controle:IncompatibleDims', 'O número de linhas de B deve ser igual ao de A (nx = %d).', nx);
end
if dims_E(1) ~= nx
    error('hinf_controle:IncompatibleDims', 'O número de linhas de E deve ser igual ao de A (nx = %d).', nx);
end
if dims_Cz(2) ~= nx
    error('hinf_controle:IncompatibleDims', 'O número de colunas de Cz deve ser igual ao de A (nx = %d).', nx);
end
if dims_Dz(1) ~= nz || dims_Dz(2) ~= nu
    error('hinf_controle:IncompatibleDims', 'Dimensões de Dz devem ser (nz x nu) -> (%d x %d).', nz, nu);
end
if dims_Ez(1) ~= nz || dims_Ez(2) ~= nw
    error('hinf_controle:IncompatibleDims', 'Dimensões de Ez devem ser (nz x nw) -> (%d x %d).', nz, nw);
end

% 4. Valida conteúdo numérico e constância das dimensões em TODOS os vértices
for k = 1:num_vertices
    if ~isnumeric(A{k}) || ~isnumeric(B{k}) || ~isnumeric(E{k}) || ...
       ~isnumeric(Cz{k}) || ~isnumeric(Dz{k}) || ~isnumeric(Ez{k})
        error('hinf_controle:NonNumeric', 'O conteúdo de todas as células no vértice %d deve ser numérico.', k);
    end
    
    if ~isequal(size(A{k}), dims_A) || ~isequal(size(B{k}), dims_B) || ...
       ~isequal(size(E{k}), dims_E) || ~isequal(size(Cz{k}), dims_Cz) || ...
       ~isequal(size(Dz{k}), dims_Dz) || ~isequal(size(Ez{k}), dims_Ez)
        error('hinf_controle:VaryingDims', 'As dimensões das matrizes não podem variar entre os vértices (verifique o vértice %d).', k);
    end
end
%   Fim do Bloco de Validação de Entrada  

%   Tratamento das Opções (varargin)  
options = [];
if nargin > 4
    if nargin == 5 && isstruct(varargin{1})
        options = varargin{1};
    else
        try
            options = struct(varargin{:});
        catch ME
            error('hinf_cost_d_rol_v2:InvalidOptions', 'Opções inválidas. Use struct ou pares nome-valor.');
        end
    end
end

% Definição de Valores Padrão
if ~isfield(options, 'solver'), options.solver = 'mosek'; end
if ~isfield(options, 'deg'), options.deg = 1; end
if ~isfield(options, 'degGamma'), options.degGamma = 0; end
if ~isfield(options, 'op'), options.op = 0; end
if ~isfield(options, 'varFolga'), options.varFolga = 1; end
if ~isfield(options, 'verbose'), options.verbose = 0; end

% Validação das Opções
if ~ischar(options.solver), error('Solver deve ser string.'); end
if options.deg < 0, error('Grau deve ser não-negativo.'); end

%   Início da Lógica da Função  

out.cpusec_m = clock;

% Preserva as matrizes originais para a validação final
Ao = A; Bo = B; Eo = E; Czo = Cz; Dzo = Dz; Ezo = Ez;

% Converte para variáveis ROLMIP (dependentes de parâmetros afins)
A_rol = rolmipvar(A, 'A', num_vertices, 1);
B_rol = rolmipvar(B, 'B', num_vertices, 1);
E_rol = rolmipvar(E, 'E', num_vertices, 1);
Cz_rol = rolmipvar(Cz, 'Cz', num_vertices, 1);
Dz_rol = rolmipvar(Dz, 'Dz', num_vertices, 1);
Ez_rol = rolmipvar(Ez, 'Ez', num_vertices, 1);

% Gera a variável de Lyapunov P(alpha)
X = rolmipvar(nx, nx, 'X', 'full', num_vertices, 0);
W = rolmipvar(nx, nx, 'W', 'symmetric', num_vertices, options.deg);
Z = rolmipvar(nu, nx, 'Z', 'full', num_vertices, 0);
csi = 0;

% Configuração da variável de desempenho (rho/gamma)
if options.op == 1
    % Caso op=1: rho é uma variável escalar constante (pico)
    rho = sdpvar(1, 1);
    rho_alpha = rolmipvar(rho, 'rho', num_vertices, 0); % Converte para rolmipvar grau 0
    obj = rho;
else        
    % Caso op=0: rho é um polinômio rho(alpha)
    powers = gen_coefs(num_vertices, options.degGamma);
    obj = [];
    rhos = cell(1, size(powers, 1));
    
    for i = 1:size(powers, 1)
        rhos{i} = sdpvar(1, 1);
        % Coeficiente binomial para ponderação na função objetivo (integral/média)
        c_val = factorial(options.degGamma) / prod(factorial(powers(i, :)));
        
        if isempty(obj)
            obj = rhos{i} * (1/c_val);
        else
            obj = rhos{i} * (1/c_val) + obj;
        end
    end    
    rho_alpha = rolmipvar(rhos, 'rho', num_vertices, options.degGamma);
end

% Montagem das LMIs    
LMIs = [W >= 0];

M11 = csi*A_rol*X + csi*B_rol*Z + csi*X'*A_rol' + csi*Z'*B_rol' - W;
M21 = X'*A_rol' + Z'*B_rol' - csi*X;
M31 = csi*Cz_rol*X + csi*Dz_rol*Z;
M41 = zeros(nw, nx);
M22 = W - X - X';
M32 = Cz_rol*X + Dz_rol*Z;
M42 = E_rol';
M33 = -rho_alpha*eye(nz);
M43 = Ez_rol';
M44 = -eye(nw);
M = [M11, M21', M31', M41';
     M21, M22, M32', M42';
     M31, M32, M33, M43';
     M41, M42, M43, M44];

LMIs = [LMIs, M <= 0];  

% Informações sobre o problema
info = lmiinfo(LMIs);
out.L = sum(info.sdp(:, 1));
out.V = size(getvariables(LMIs), 2);

% Tempo de montagem
out.cpusec_m = etime(clock, out.cpusec_m);

% Solução das LMIs
sol = solvesdp(LMIs, obj, sdpsettings('verbose', options.verbose, 'solver', options.solver));

% Tempo de solução
out.cpusec_s = sol.solvertime;

% Verificação de viabilidade (residual primal mínimo)
delta = min(checkset(LMIs));

out.feas = 0;
out.wc = Inf;
out.normSysWc = Inf;
out.alpha = [];
out.gcosts = [];
out.realCosts = [];
out.wc_alpha = [];
out.Z_cells = [];
out.X_cells = []; 
out.K_wc = [];

if delta > -1e-5
    out.feas = 1;
    
    % Extrai os valores numéricos (células dos vértices) de Z e X
    out.Z_cells = double(Z);
    out.X_cells = double(X); 
    
    %   Validação (Grid Search em todo o Simplex) 
    maxNorm = -1e10;
    passo_grid = 0.05; 
    
    % Gera os pontos do simplex unitário usando a nova função
    % NOTA: Se N=1, a função retorna {[1]} e o loop roda perfeitamente 1 vez!
    pontos_alpha = particao_simplex(num_vertices, passo_grid);
    out.alpha = zeros(length(pontos_alpha), num_vertices);
    
    for k = 1:length(pontos_alpha)
        alpha_vec = pontos_alpha{k}; % Vetor coluna
        s = alpha_vec'; % Vetor linha para o ROLMIP
        out.alpha(k, :) = s;
        
        % 1. Custo Garantido (via LMI)
        val = double(evalpar(rho_alpha, {s}));
        v = sqrt(max(0, val)); % max(0,.) evita complexos por ruído numérico
        out.gcosts = [out.gcosts; v];
        
        % Rastreia o pior caso garantido
        if v > maxNorm
            maxNorm = v;
            wc_vec = s;
        end
        
        % 2. Custo Real (Sistema em MALHA FECHADA Congelado)
        Z_val = double(evalpar(Z, {s}));
        X_val = double(evalpar(X, {s}));
        K_val = Z_val / X_val; 
        
        A_val  = double(evalpar(A_rol, {s}));
        B_val  = double(evalpar(B_rol, {s}));
        E_val  = double(evalpar(E_rol, {s}));
        Cz_val = double(evalpar(Cz_rol, {s}));
        Dz_val = double(evalpar(Dz_rol, {s}));
        Ez_val = double(evalpar(Ez_rol, {s}));
        
        Acl_val = A_val + B_val * K_val;
        Bcl_val = E_val;
        Ccl_val = Cz_val + Dz_val * K_val;
        Dcl_val = Ez_val;
        
        sys_local = ss(Acl_val, Bcl_val, Ccl_val, Dcl_val, -1);
        try
            out.realCosts = [out.realCosts; norm(sys_local, inf)];               
        catch
            out.realCosts = [out.realCosts; NaN];
        end
    end
    
    out.wc_alpha = wc_vec;
    out.wc = maxNorm; 
    
    % Recalcula o sistema e o controlador no ponto de PIOR CASO encontrado
    Z_wc = double(evalpar(Z, {wc_vec}));
    X_wc = double(evalpar(X, {wc_vec}));
    
    out.K_wc = Z_wc / X_wc;
    
    % Salva as matrizes de malha fechada do pior caso na estrutura de saída
    out.Acl_wc = double(evalpar(A_rol, {wc_vec})) + double(evalpar(B_rol, {wc_vec})) * out.K_wc;
    out.Bcl_wc = double(evalpar(E_rol, {wc_vec}));
    out.Ccl_wc = double(evalpar(Cz_rol, {wc_vec})) + double(evalpar(Dz_rol, {wc_vec})) * out.K_wc;
    out.Dcl_wc = double(evalpar(Ez_rol, {wc_vec}));
    
    % Calcula a norma real do pior caso
    sys_wc = ss(out.Acl_wc, out.Bcl_wc, out.Ccl_wc, out.Dcl_wc, -1);
    try
        out.normSysWc = norm(sys_wc, inf);
    catch
        out.normSysWc = NaN;
    end
end

end