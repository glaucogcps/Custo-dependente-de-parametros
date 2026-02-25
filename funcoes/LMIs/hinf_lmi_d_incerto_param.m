function out = hinf_lmi_d_incerto_param(A, B, C, D, varargin)
% function out = hinf_lmi_d_incerto_param(A, B, C, D, varargin)
%
% Avalia o custo H infinito garantido de um sistema linear discreto no tempo 
% com incerteza politópica, utilizando funções de Lyapunov dependentes de parâmetros.
%
% Esta função implementa duas abordagens para o cálculo do custo:
% 1. Lema Real Limitado (Bounded Real Lemma) clássico dependente de parâmetros.
% 2. Uma formulação relaxada utilizando variáveis de folga (multiplicadores),
%    que desacopla as matrizes do sistema da matriz de Lyapunov (Lema de Finsler/Projeção).
%
% input:
%   A, B, C, D -> Cell arrays contendo as matrizes de espaço de estados (vértices do politopo).
%                 Cada cell array deve conter N matrizes.
%   varargin   -> (opcional) Pares nome-valor ou estrutura de opções:
%                 - 'solver': String. Nome do solver (padrão: 'sedumi').
%                 - 'deg': Inteiro. Grau do polinômio da matriz de Lyapunov P(alpha). (Padrão: 1)
%                 - 'degGamma': Inteiro. Grau do polinômio do custo gamma(alpha). (Padrão: 0)
%                 - 'op': Inteiro (0 ou 1). 
%                         0 -> Otimiza a integral/média do custo (gamma polinomial).
%                         1 -> Otimiza o pico do custo (gamma constante). (Padrão: 0)
%                 - 'varFolga': Booleano (0 ou 1).
%                         0 -> Formulação clássica (BRL).
%                         1 -> Formulação relaxada com variáveis de folga. (Padrão: 1)
%                 - 'verbose': Inteiro. Nível de verbosidade do solver. (Padrão: 0)
%
% output:
%   out.feas       -> 1 se viável, 0 se inviável.
%   out.wc         -> Pior caso (worst-case) da norma Hinf encontrado na validação.
%   out.normSysWc  -> Norma Hinf real do sistema congelado no ponto de pior caso.
%   out.cpusec_m   -> Tempo de montagem das LMIs.
%   out.cpusec_s   -> Tempo de solução do solver.
%   out.L          -> Número de linhas de LMI.
%   out.V          -> Número de variáveis escalares.
%   out.gcosts     -> Vetor de custos garantidos avaliados no grid de validação.
%   out.realCosts  -> Vetor de custos reais avaliados no grid de validação.
%
% Date: 28/01/2026

%   Bloco de Validação de Entrada  
if nargin < 4
    error('hinf_cost_d_rol_v2:NotEnoughInputs', 'São necessários A, B, C, D (cell arrays de matrizes).');
end

input_cells = {A, B, C, D};
input_names = {'A', 'B', 'C', 'D'};
num_vertices = numel(A); % N

% Valida tipos e consistência do número de vértices
for i = 1:length(input_cells)
    if ~iscell(input_cells{i})
        error('hinf_cost_d_rol_v2:InvalidInputType', 'O argumento %s deve ser um cell array.', input_names{i});
    end
    if numel(input_cells{i}) ~= num_vertices
        error('hinf_cost_d_rol_v2:DimensionMismatch', 'Todos os cell arrays devem ter o mesmo número de vértices.');
    end
end

% Valida dimensões das matrizes
n = size(A{1}, 1); % Ordem do sistema
m = size(B{1}, 2); % Entradas
p = size(C{1}, 1); % Saídas

dims_A = size(A{1});
dims_B = size(B{1});
dims_C = size(C{1});
dims_D = size(D{1});

if dims_A(1) ~= dims_A(2)
    error('hinf_cost_d_rol_v2:InvalidADim', 'Matrizes A devem ser quadradas.');
end
if dims_B(1) ~= n || dims_C(2) ~= n
    error('hinf_cost_d_rol_v2:IncompatibleDims', 'Dimensões de B ou C incompatíveis com A.');
end
if dims_D(1) ~= p || dims_D(2) ~= m
    error('hinf_cost_d_rol_v2:IncompatibleDims', 'Dimensões de D incompatíveis com C e B.');
end

% Valida conteúdo numérico
for k = 1:num_vertices
    if ~isnumeric(A{k}) || ~isnumeric(B{k}) || ~isnumeric(C{k}) || ~isnumeric(D{k})
        error('hinf_cost_d_rol_v2:NonNumeric', 'Conteúdo das células deve ser numérico.');
    end
    if ~isequal(size(A{k}), dims_A) || ~isequal(size(B{k}), dims_B) || ...
       ~isequal(size(C{k}), dims_C) || ~isequal(size(D{k}), dims_D)
        error('hinf_cost_d_rol_v2:VaryingDims', 'Dimensões das matrizes variam entre os vértices.');
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
    %   Formulação Clássica (Bounded Real Lemma)  
    % [ A'P A + P A + C'C    P B + C'D      ]
    % [ B'P + D'C            D'D - rho*I    ] < 0
    % Nota: O código original usava uma estrutura discreta/contínua mista ou específica.
    % Abaixo segue a lógica estrita do código original fornecido:
    
    T11 = A_rol' * P + P * A_rol + C_rol' * C_rol; 
    
    T21 = B_rol' * P + D_rol' * C_rol;
    T22 = D_rol' * D_rol - rho_alpha * eye(m);
    
    T = [T11, T21'; T21, T22];
    
    LMIs = [T <= 0, P >= 0];
    
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
out.wc = [];
out.normSysWc = [];
out.alpha = [];
out.gcosts = [];
out.realCosts = [];
out.wc_alpha = [];

if delta > -1e-5
    out.feas = 1;
    
    %   Validação (Grid Search nas Arestas)  
    % Varre as arestas entre vértices para comparar o custo garantido vs. real
    maxNorm = -1e10;
    
    for i = 1:num_vertices
        for j = i+1:num_vertices
            for a = 0:0.01:1
                % Vetor de coordenadas baricêntricas para a aresta i-j
                s = zeros(1, num_vertices);
                s(1, i) = a;
                s(1, j) = 1 - a;
                
                out.alpha = [out.alpha; a];
                
                % 1. Custo Garantido (via LMI)
                v = double(sqrt(evalpar(rho_alpha, {s})));
                out.gcosts = [out.gcosts; v];
                
                % Rastreia o pior caso garantido
                if v > maxNorm
                    maxNorm = v;
                    wc_vec = s;
                end
                
                % 2. Custo Real (Sistema Congelado)
                Aa_val = a * Ao{i} + (1 - a) * Ao{j};
                Ba_val = a * Bo{i} + (1 - a) * Bo{j};
                Ca_val = a * Co{i} + (1 - a) * Co{j};
                Da_val = a * Do{i} + (1 - a) * Do{j};
                
                % Cria sistema SS (assumindo tempo de amostragem -1 para discreto genérico)
                sys_local = ss(Aa_val, Ba_val, Ca_val, Da_val, -1);
                out.realCosts = [out.realCosts; norm(sys_local, inf)];                
            end
        end
    end
    
    out.wc_alpha = wc_vec;
    out.wc = maxNorm; 
    
    % Recalcula o sistema no ponto de pior caso encontrado
    Aa_wc = zeros(n); Ba_wc = zeros(n, m); Ca_wc = zeros(p, n); Da_wc = zeros(p, m);
    for i = 1:num_vertices
        Aa_wc = Aa_wc + wc_vec(1, i) * Ao{i};
        Ba_wc = Ba_wc + wc_vec(1, i) * Bo{i};
        Ca_wc = Ca_wc + wc_vec(1, i) * Co{i};
        Da_wc = Da_wc + wc_vec(1, i) * Do{i};
    end
    
    sys_wc = ss(Aa_wc, Ba_wc, Ca_wc, Da_wc, -1);
    out.normSysWc = norm(sys_wc, inf);
end

end