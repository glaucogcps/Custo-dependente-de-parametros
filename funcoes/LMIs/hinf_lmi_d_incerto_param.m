function out = hinf_lmi_d_incerto_param(A, B, C, D, varargin)
% function out = hinf_lmi_d_incerto_param(A, B, C, D, varargin)
%
% Avalia o custo H infinito garantido de um SISTEMA LINEAR DISCRETO no tempo 
% com incerteza politópica, utilizando funções de Lyapunov dependentes de parâmetros.
%
% Esta função implementa duas abordagens para o cálculo do custo:
% 1. Lema Real Limitado (Bounded Real Lemma) clássico dependente de parâmetros.
% 2. Formulação relaxada utilizando variáveis de folga (multiplicadores),
%    que desacopla as matrizes do sistema da matriz de Lyapunov.
%
% Entradas:
%   A, B, C, D - Cell arrays contendo as matrizes dos vértices do politopo.
%                Ex: A = {A1, A2, ..., AN}, onde N é o número de vértices.
%   varargin   - (Opcional) Struct ou pares nome-valor com opções:
%                .solver   : Solver a ser utilizado (default: 'mosek')
%                .deg      : Grau do polinômio da matriz de Lyapunov P (default: 1)
%                .degGamma : Grau do polinômio da variável de desempenho (default: 0)
%                .op       : Tipo de otimização (0 = integral/média, 1 = pico) (default: 0)
%                .varFolga : Uso de variáveis de folga/Finsler (1 = sim, 0 = não) (default: 1)
%                .verbose  : Nível de detalhe do solver (default: 0)
%
% Saída:
%   out - Struct contendo os resultados da otimização e validação:
%         .feas      : Status de viabilidade (1 = viável, 0 = inviável)
%         .wc        : Pior custo Hinf garantido encontrado no grid
%         .normSysWc : Custo Hinf real do sistema no ponto de pior caso garantido
%         .wc_alpha  : Coordenadas baricêntricas do pior caso
%         .gcosts    : Vetor com os custos garantidos avaliados no grid
%         .realCosts : Vetor com os custos reais (norma inf) avaliados no grid
%         .alpha     : Matriz com os pontos do simplex testados (cada linha é um ponto)
%         .L         : Número de linhas das LMIs
%         .V         : Número de variáveis escalares
%         .cpusec_m  : Tempo de montagem das LMIs
%         .cpusec_s  : Tempo de resolução do solver
%         .rho_alpha : Objeto ROLMIP do polinômio de custo (gamma^2)
%
% Date: 22/03/2026

%   Bloco de Validação de Entrada  
if nargin < 4
    error('hinf_lmi_d_param:NotEnoughInputs', 'São necessários A, B, C, D (cell arrays de matrizes).');
end

input_cells = {A, B, C, D};
input_names = {'A', 'B', 'C', 'D'};
num_vertices = numel(A); % N

% Valida tipos e consistência do número de vértices
for i = 1:length(input_cells)
    if ~iscell(input_cells{i})
        error('hinf_lmi_d_param:InvalidInputType', 'O argumento %s deve ser um cell array.', input_names{i});
    end
    if numel(input_cells{i}) ~= num_vertices
        error('hinf_lmi_d_param:DimensionMismatch', 'Todos os cell arrays devem ter o mesmo número de vértices.');
    end
end

% Valida dimensões das matrizes
n = size(A{1}, 1); % Ordem do sistema
m = size(B{1}, 2); % Entradas
p = size(C{1}, 1); % Saídas

dims_A = size(A{1}); dims_B = size(B{1}); dims_C = size(C{1}); dims_D = size(D{1});

if dims_A(1) ~= dims_A(2)
    error('hinf_lmi_d_param:InvalidADim', 'Matrizes A devem ser quadradas.');
end
if dims_B(1) ~= n || dims_C(2) ~= n
    error('hinf_lmi_d_param:IncompatibleDims', 'Dimensões de B ou C incompatíveis com A.');
end
if dims_D(1) ~= p || dims_D(2) ~= m
    error('hinf_lmi_d_param:IncompatibleDims', 'Dimensões de D incompatíveis com C e B.');
end

% Valida conteúdo numérico
for k = 1:num_vertices
    if ~isnumeric(A{k}) || ~isnumeric(B{k}) || ~isnumeric(C{k}) || ~isnumeric(D{k})
        error('hinf_lmi_d_param:NonNumeric', 'Conteúdo das células deve ser numérico.');
    end
    if ~isequal(size(A{k}), dims_A) || ~isequal(size(B{k}), dims_B) || ...
       ~isequal(size(C{k}), dims_C) || ~isequal(size(D{k}), dims_D)
        error('hinf_lmi_d_param:VaryingDims', 'Dimensões das matrizes variam entre os vértices.');
    end
end

%   Tratamento das Opções (varargin)  
options = struct('solver', 'mosek', 'deg', 1, 'degGamma', 0, 'op', 0, 'varFolga', 1, 'verbose', 0);
if nargin > 4
    if nargin == 5 && isstruct(varargin{1})
        opts_in = varargin{1};
    else
        try opts_in = struct(varargin{:}); catch, opts_in = struct(); end
    end
    f = fieldnames(opts_in);
    for i = 1:length(f), options.(f{i}) = opts_in.(f{i}); end
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
Ao = A; Bo = B; Co = C; Do = D;

% Converte para variáveis ROLMIP (dependentes de parâmetros afins)
A_rol = rolmipvar(A, 'A', num_vertices, 1);
B_rol = rolmipvar(B, 'B', num_vertices, 1);
C_rol = rolmipvar(C, 'C', num_vertices, 1);
D_rol = rolmipvar(D, 'D', num_vertices, 1);

% Gera a variável de Lyapunov P(alpha)
P = rolmipvar(n, n, 'P', 'symmetric', num_vertices, options.deg);

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
if options.varFolga == 0
    M11 = A_rol'*P*A_rol - P + C_rol'*C_rol;
    M21 = B_rol'*P*A_rol + D_rol'*C_rol;
    M22 = D_rol'*D_rol + B_rol'*P*B_rol - rho_alpha * eye(m);
    M = [M11, M21';
         M21, M22];
    LMIs = [M <= 0, P >= 0];
    
else
    %   Formulação Relaxada (Variáveis de Folga / Finsler)  
    % Introduz multiplicador X para desacoplar as matrizes do sistema de P
    
    X = rolmipvar(2*n + p + m, n + p, 'X', 'full', num_vertices, options.deg);
    Z = zeros(n);
    
    % Matriz Q contendo P e rho
    Q = blkdiag([P, Z; Z, -P], eye(p), -rho_alpha * eye(m));
    
    % Matriz Bcal contendo a dinâmica do sistema
    % Estrutura típica para impor a dinâmica via multiplicadores
    Bcal = [eye(n), -A_rol, zeros(n, p), -B_rol;
            zeros(p, n), -C_rol, eye(p), -D_rol];
            
    % LMI: Q + X*Bcal + Bcal'*X' <= 0
    LMIs = [Q + X * Bcal + Bcal' * X' <= 0, P >= 0];
end

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
out.rho_alpha = [];

if delta > -1e-5
    out.feas = 1;
    out.rho_alpha = rho_alpha;
    
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
        val = double(evalpar(rho_alpha, {s}));
        v = sqrt(max(0, val)); % Raiz quadrada pois otimizamos gamma^2
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
        
        % Cria sistema SS (assumindo tempo de amostragem -1 para discreto genérico)
        sys_local = ss(Aa_val, Ba_val, Ca_val, Da_val, -1);
        try
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
    
    sys_wc = ss(Aa_wc, Ba_wc, Ca_wc, Da_wc, -1);
    try
        out.normSysWc = norm(sys_wc, inf);
    catch
        out.normSysWc = NaN;
    end
end

end