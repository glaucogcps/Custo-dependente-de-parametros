function out = hinf_lmi_c_incerto_param(A, B, C, D, varargin)
% function out = hinf_lmi_c_incerto_param(A, B, C, D, varargin)
%
% Avalia o custo H infinito garantido de um SISTEMA LINEAR CONTÍNUO
% com incerteza politópica, utilizando funções de Lyapunov dependentes de parâmetros.
%
% Entradas:
%   A, B, C, D - Cell arrays contendo as matrizes dos vértices do politopo.
%                Ex: A = {A1, A2, ..., AN}, onde N é o número de vértices.
%   varargin   - (Opcional) Struct com opções de configuração:
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
    error('hinf_lmi_c_param:NotEnoughInputs', 'São necessários A, B, C, D (cell arrays de matrizes).');
end

input_cells = {A, B, C, D};
input_names = {'A', 'B', 'C', 'D'};
num_vertices = numel(A); % N

% Valida tipos e consistência do número de vértices
for i = 1:length(input_cells)
    if ~iscell(input_cells{i})
        error('hinf_lmi_c_param:InvalidInputType', 'O argumento %s deve ser um cell array.', input_names{i});
    end
    if numel(input_cells{i}) ~= num_vertices
        error('hinf_lmi_c_param:DimensionMismatch', 'Todos os cell arrays devem ter o mesmo número de vértices.');
    end
end

% Valida dimensões das matrizes
n = size(A{1}, 1); m = size(B{1}, 2); p = size(C{1}, 1);
dims_A = size(A{1}); dims_B = size(B{1}); dims_C = size(C{1}); dims_D = size(D{1});

if dims_A(1) ~= dims_A(2)
    error('hinf_lmi_c_param:InvalidADim', 'Matrizes A devem ser quadradas.');
end
if dims_B(1) ~= n || dims_C(2) ~= n
    error('hinf_lmi_c_param:IncompatibleDims', 'Dimensões de B ou C incompatíveis com A.');
end
if dims_D(1) ~= p || dims_D(2) ~= m
    error('hinf_lmi_c_param:IncompatibleDims', 'Dimensões de D incompatíveis com C e B.');
end

% Valida conteúdo numérico
for k = 1:num_vertices
    if ~isnumeric(A{k}) || ~isnumeric(B{k}) || ~isnumeric(C{k}) || ~isnumeric(D{k})
        error('hinf_lmi_c_param:NonNumeric', 'Conteúdo das células deve ser numérico.');
    end
    if ~isequal(size(A{k}), dims_A) || ~isequal(size(B{k}), dims_B) || ...
       ~isequal(size(C{k}), dims_C) || ~isequal(size(D{k}), dims_D)
        error('hinf_lmi_c_param:VaryingDims', 'Dimensões das matrizes variam entre os vértices.');
    end
end

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

%   Variável de Desempenho   
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

if options.varFolga == 0
    M11 = A_rol'*P + P*A_rol + C_rol'*C_rol;
    M21 = B_rol'*P + D_rol'*C_rol;
    M22 = D_rol'*D_rol - rho_alpha * eye(m);
    M = [M11, M21';
         M21, M22];
    LMIs = [M <= 0, P >= 0];
    
else
    %   Formulação Relaxada (Variáveis de Folga / Finsler)  
    % Introduz multiplicador X para desacoplar as matrizes do sistema de P
    
    X = rolmipvar(2*n + p + m, n + p, 'X', 'full', num_vertices, options.deg);
    Z = zeros(n);
    
    % Matriz Q contendo P e rho
    Q = blkdiag([Z, P; P, Z], eye(p), -rho_alpha * eye(m));
    
    % Matriz Bcal contendo a dinâmica do sistema
    % Estrutura típica para impor a dinâmica via multiplicadores
    Bcal = [eye(n), -A_rol, zeros(n, p), -B_rol;
            zeros(p, n), -C_rol, eye(p), -D_rol];
            
    % LMI: Q + X*Bcal + Bcal'*X' <= 0
    LMIs = [Q + X * Bcal + Bcal' * X' <= 0, P >= 0];
end

%   Estatísticas de Complexidade  
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
    
    % Inicializa vetores para plotagem
    max_val = -inf;
    out.gcosts = [];
    out.realCosts = [];
    
    % Define o passo da malha (pode ser parametrizado no varargin futuramente)
    passo_grid = 0.05; 
    
    % Gera os pontos do simplex unitário usando a nova função
    pontos_alpha = particao_simplex(num_vertices, passo_grid);
    
    % Inicializa a matriz para salvar os pontos testados
    out.alpha = zeros(length(pontos_alpha), num_vertices);
    
    % Preserva originais para validação
    Ao = A; Bo = B; Co = C; Do = D;
    
    for k = 1:length(pontos_alpha)
        alpha_vec = pontos_alpha{k}; % Vetor coluna
        s = alpha_vec'; % Vetor linha para o ROLMIP
        out.alpha(k, :) = s;
        
        % 1. Custo Garantido H-infinito (rho_alpha representa gamma^2)
        val = double(evalpar(rho_alpha, {s}));
        v = sqrt(max(0, val)); % Raiz quadrada pois otimizamos gamma^2
        out.gcosts = [out.gcosts; v];
        
        if v > max_val, max_val = v; end
        
        % 2. Custo Real H-infinito
        Aa = zeros(n, n); Ba = zeros(n, m);
        Ca = zeros(p, n); Da = zeros(p, m);
        
        % Combinação convexa generalizada para N vértices
        for i = 1:num_vertices
            Aa = Aa + alpha_vec(i) * Ao{i};
            Ba = Ba + alpha_vec(i) * Bo{i};
            Ca = Ca + alpha_vec(i) * Co{i};
            Da = Da + alpha_vec(i) * Do{i};
        end
        
        % Cria o sistema contínuo
        sys = ss(Aa, Ba, Ca, Da);
        
        try
            % Calcula a norma H-infinito real do sistema
            out.realCosts = [out.realCosts; norm(sys, inf)];
        catch
            out.realCosts = [out.realCosts; NaN];
        end
    end
    out.wc = max_val; 
else
    out.feas = 0;
    out.rho_alpha = [];
    out.wc = Inf;
    out.gcosts = [];
    out.realCosts = [];
    out.alpha = [];
end
end