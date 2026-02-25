function out = h2_lmi_d_incerto(A,B,C,D, varargin)
% function out = h2_lmi_d_incerto(A,B,C,D, varargin)
%
% Avalia a norma H2 de um sistema linear discreto no tempo com incerteza politópica.
% As Desigualdades Matriciais Lineares (LMIs) são programadas usando YALMIP e ROLMIP,
% e podem ser resolvidas por qualquer solver LMI suportado pelo YALMIP.
%
% Esta função busca a menor norma H2 para o sistema incerto, considerando
% uma representação politópica para as matrizes de espaço de estados e
% utilizando funções de Lyapunov polinomiais (se 'grau' > 0).
%
% input:
%   A, B, C, D -> Cell arrays contendo as matrizes de espaço de estados (vértices do politopo).
%                 Cada cell array deve conter N matrizes, onde N é o número de vértices do politopo.
%                 As dimensões das matrizes devem ser compatíveis para um sistema de espaço de estados:
%                 A: n x n, B: n x m, C: p x n, D: p x m.
%   varargin   -> (opcional) Estrutura de opções ou pares nome-valor para configurar a rotina.
%                 Campos suportados:
%                 - 'solver': String. Nome do solver LMI a ser utilizado (e.g., 'mosek', 'sedumi').
%                             (Padrão: 'mosek')
%                 - 'h2': Scalar numérico. Se diferente de 0, a rotina tenta verificar a
%                         viabilidade de uma norma H2 menor ou igual a este valor.
%                         Se 0, a rotina otimiza para encontrar a menor norma H2.
%                         (Padrão: 0)
%                 - 'tolerance': Scalar numérico. Tolerância para a verificação de viabilidade das LMIs.
%                                (Padrão: 1e-7)
%                 - 'grau': Scalar inteiro. Grau do polinômio para as matrizes de Lyapunov (P e M).
%                           'grau = 0' implica matrizes constantes.
%                           'grau = 1' implica matrizes afins (lineares em alpha).
%                           (Padrão: 1)
%
% output:
%   out.h2       -> Norma H2 calculada (0 se o problema for infactível ou se 'h2' for 0 e a otimização falhar).
%   out.cpusec   -> Tempo de CPU para resolver as LMIs (segundos).
%   out.cpusec_m -> Tempo de CPU para montar as LMIs (segundos).
%   out.P        -> Matriz de Lyapunov P (polinomial, avaliada no ponto ótimo).
%   out.M        -> Matriz Gramiana M (polinomial, avaliada no ponto ótimo).
%   out.V        -> Número de variáveis de decisão no problema YALMIP.
%   out.L        -> Número de linhas de LMIs no problema YALMIP.
%   out.p        -> Valor mínimo do residual das LMIs (resultado de checkset), indicando a viabilidade.
%
% Date: 25/09/2025
% Author: glauco.gcps@gmail.com

%   Bloco de Validação de Entrada  
if nargin < 4
    error('h2_lmi_d_incerto:NotEnoughInputs', 'São necessários A, B, C, D (cell arrays de matrizes).');
end

input_cells = {A, B, C, D};
input_names = {'A', 'B', 'C', 'D'};
num_vertices = numel(A); % N, o número de vértices do politopo

% Valida que todas as entradas são cell arrays e têm o mesmo número de elementos
for i = 1:length(input_cells)
    if ~iscell(input_cells{i})
        error('h2_lmi_d_incerto:InvalidInputType', 'O argumento %s deve ser um cell array.', input_names{i});
    end
    if numel(input_cells{i}) ~= num_vertices
        error('h2_lmi_d_incerto:DimensionMismatch', 'Todos os cell arrays (A, B, C, D) devem ter o mesmo número de vértices (elementos).');
    end
end

% Valida as dimensões e tipos das matrizes dentro dos cell arrays
order = size(A{1},1); % n, ordem do sistema
inputs = size(B{1},2); % m, número de entradas
outputs = size(C{1},1); % p, número de saídas

first_A_dims = size(A{1});
first_B_dims = size(B{1});
first_C_dims = size(C{1});
first_D_dims = size(D{1});

% Verifica a compatibilidade para o sistema de espaço de estados (A,B,C,D)
% A: n x n, B: n x m, C: p x n, D: p x m
if first_A_dims(1) ~= first_A_dims(2)
    error('h2_lmi_d_incerto:InvalidADim', 'Matrizes A devem ser quadradas (n x n).');
end
if first_B_dims(1) ~= first_A_dims(1)
    error('h2_lmi_d_incerto:IncompatibleDims', 'Dimensões de B (linhas) incompatíveis com A (linhas).');
end
if first_C_dims(2) ~= first_A_dims(1)
    error('h2_lmi_d_incerto:IncompatibleDims', 'Dimensões de C (colunas) incompatíveis com A (colunas).');
end
if first_D_dims(1) ~= first_C_dims(1) || first_D_dims(2) ~= first_B_dims(2)
    error('h2_lmi_d_incerto:IncompatibleDims', 'Dimensões de D incompatíveis com C e B.');
end

for k = 1:num_vertices
    current_A = A{k};
    current_B = B{k};
    current_C = C{k};
    current_D = D{k};

    if ~ismatrix(current_A) || ~isnumeric(current_A) || ~ismatrix(current_B) || ~isnumeric(current_B) || ...
       ~ismatrix(current_C) || ~isnumeric(current_C) || ~ismatrix(current_D) || ~isnumeric(current_D)
        error('h2_lmi_d_incerto:InvalidMatrixElement', 'Todos os elementos dos cell arrays devem ser matrizes numéricas.');
    end
    
    if ~isequal(size(current_A), first_A_dims) || ~isequal(size(current_B), first_B_dims) || ...
       ~isequal(size(current_C), first_C_dims) || ~isequal(size(current_D), first_D_dims)
        error('h2_lmi_d_incerto:DimensionMismatch', 'Todas as matrizes dentro de cada cell array devem ter as mesmas dimensões.');
    end
end
%   Fim do Bloco de Validação de Entrada  

%   Tratamento das Opções (varargin)  
options = [];
if nargin > 4 % A, B, C, D são 4 argumentos obrigatórios
    if nargin == 5 && isstruct(varargin{1}) % A, B, C, D, options_struct
        options = varargin{1};
    else % A, B, C, D, 'option_name', option_value, ...
        try
            options = struct(varargin{:});
        catch ME
            error('h2_lmi_d_incerto:InvalidOptions', 'Argumentos opcionais inválidos. Use uma estrutura ou pares nome-valor. Erro: %s', ME.message);
        end
    end
end

%   Definição de Valores Padrão para as Opções  
if ~isfield(options,'solver')
    options.solver = 'mosek';
end
if ~isfield(options,'h2')
    options.h2 = 0;
end
if ~isfield(options,'tolerance')
    options.tolerance = 1e-7;
end
if ~isfield(options,'grau')
    options.grau = 1; % Padrão: Matrizes de Lyapunov afins (grau 1)
end

grau_mu = 0;
if isfield(options,'grau_mu')
    grau_mu = options.grau_mu; 
end

%   Validação das Opções  
if ~ischar(options.solver) || isempty(options.solver)
    error('h2_lmi_d_incerto:InvalidSolver', 'A opção ''solver'' deve ser uma string não vazia.');
end
if ~isnumeric(options.h2) || ~isscalar(options.h2) || options.h2 < 0
    error('h2_lmi_d_incerto:InvalidH2Value', 'A opção ''h2'' deve ser um escalar numérico não negativo.');
end
if ~isnumeric(options.tolerance) || ~isscalar(options.tolerance) || options.tolerance <= 0
    error('h2_lmi_d_incerto:InvalidTolerance', 'A opção ''tolerance'' deve ser um escalar numérico positivo.');
end
if ~isnumeric(options.grau) || ~isscalar(options.grau) || options.grau < 0 || mod(options.grau, 1) ~= 0
    error('h2_lmi_d_incerto:InvalidGrau', 'A opção ''grau'' deve ser um inteiro não negativo.');
end

grau = options.grau;

%   Início da Lógica da Função  
% As variáveis 'vertices', 'order', 'inputs', 'outputs' já foram determinadas
% no bloco de validação de entrada.

out.cpusec_m = clock;

% Criação das variáveis YALMIP
if options.h2 == 0
    if grau_mu == 0
        mu = sdpvar(1); % Variável de otimização para a norma h2^2
        % 1. Crie a matriz numérica ou sdpvar para o termo constante.
        mu_matrix_term = -mu * eye(inputs); 
        % 2. Converta-a em um rolmipvar constante (grau 0 em alpha).
        obj = mu;
    else 
        mu = rolmipvar(1, 1, 'mu', 'full', num_vertices, grau_mu);
        soma = 0;
        m = coeffs(mu);
        for i = 1:size(m, 2)
            soma = soma + m{i};
        end
        obj = soma;
    end
else
    mu = options.h2 * options.h2; % Norma h2^2 alvo para verificação de viabilidade
    obj = []; % Sem objetivo de otimização, apenas viabilidade
    options.tolerance = 0; % Desativa a tolerância para o checkset se for apenas viabilidade
end

% Definição das matrizes do sistema como variáveis ROLMIP para incerteza politópica
% 'rolmipvar' cria variáveis que dependem do vetor de incerteza alpha
Aa = rolmipvar(A, 'Aa', num_vertices, 1); % Matriz A incerta
Ba = rolmipvar(B, 'Ba', num_vertices, 1); % Matriz B incerta
Ca = rolmipvar(C, 'Ca', num_vertices, 1); % Matriz C incerta
Da = rolmipvar(D, 'Da', num_vertices, 1); % Matriz D incerta

% Definição das matrizes de Lyapunov e Gramiana como variáveis ROLMIP
% P e M são polinomiais em alpha, com grau especificado
P = rolmipvar(order, order, 'P', 'symmetric', num_vertices, grau);
M = rolmipvar(inputs, inputs, 'M', 'symmetric', num_vertices, grau);

% LMI 1: P >= 0 e M >= 0
% Garante que as matrizes de Lyapunov e Gramiana sejam positivas definidas
LMIs = [P >= 0, M >= 0];
    
% LMI 2: Condição de traço e parte da condição de norma H2
% [M - D'D, B'P; P'B, P] >= 0 e trace(M) <= mu
% Esta LMI é uma condição para a norma H2 de sistemas discretos.
% trace(M) representa a norma H2^2 do sistema.
Q11 = M - Da'*Da;
Q12 = Ba'*P;
Q22 = P;
Q = [Q11 Q12;
    Q12' Q22];

LMIs = [LMIs, Q >= 0, trace(M) <= mu];

% LMI 3: Condição de Gramiana (Lyapunov) para estabilidade e norma H2
% [P - A'PA, C'; C, I] >= 0
% Esta LMI é equivalente a P - A'PA >= C'C, que é uma condição de estabilidade
% e parte da formulação da norma H2 para sistemas discretos.
G11 = P - Aa'*P*Aa;
G12 = Ca';
G22 = eye(outputs); % Matriz identidade de dimensão p x p
G = [G11 G12;
    G12' G22];

LMIs = [LMIs, G >= 0];

out.cpusec_m = etime(clock,out.cpusec_m);

% Determina o número de linhas de LMIs e variáveis escalares
out.L = 0;
for i=1:size(LMIs,1)
    out.L = out.L + size(sdpvar(LMIs(i)),1);
end
out.V = size(getvariables(LMIs),2);

% Resolve o problema LMI
sol = solvesdp(LMIs,obj,sdpsettings('verbose',0,'solver',options.solver));
out.cpusec = sol.solvertime;
out.p = min(checkset(LMIs)); % Verifica a viabilidade das LMIs

out.h2 = 0;
out.feas = 0;
% Captura as soluções (se existirem)
if out.p > -options.tolerance
    out.feas = 1;
    if grau_mu == 0
        out.h2 = sqrt(double(mu)); 
    else
        out.h2 = double(mu);
    end
    out.P = double(P); % Avalia a matriz P no ponto ótimo
    out.M = double(M);
end