function output = norma_h2_sistema_incerto_discrete(A_cell, B_cell, C_cell, D_cell, varargin)
% function output = norma_h2_sistema_incerto_discrete(A_cell, B_cell, C_cell, D_cell, varargin)
%
% Calcula a norma H2 de um sistema incerto discreto no tempo (A,B,C,D)(alpha).
%
% Esta função varre o simplex unitário para combinações de alpha, constrói
% o sistema (A,B,C,D)(alpha) e calcula a norma H2 para cada combinação.
% O pior caso (maior norma H2) é identificado.
%
% input:
%   A_cell, B_cell, C_cell, D_cell -> Cell arrays contendo os vértices das
%                                     matrizes A, B, C, D, respectivamente, que definem
%                                     o politopo de sistemas.
%                                     Todos os cell arrays devem ter o mesmo número de elementos (N),
%                                     representando o número de vértices do politopo.
%                                     As dimensões das matrizes devem ser compatíveis para
%                                     formar um sistema de espaço de estados (A: nxn, B: nxm, C: pxn, D: pxm).
%   varargin                       -> (opcional) Estrutura de opções ou pares nome-valor para configurar a rotina.
%                                     Campos suportados:
%                                     - 'simplex_step': Scalar numérico. Define o passo para a geração de pontos
%                                                       no simplex unitário. Um valor menor resulta em uma varredura
%                                                       mais fina e precisa do politopo, mas aumenta o tempo de computação.
%                                                       Deve ser um valor positivo e menor ou igual a 1.
%                                                       (Padrão: 0.001)
%                                     - 'Ts': Scalar numérico. Período de amostragem do sistema discreto.
%                                             Essencial para a função 'norm' interpretar corretamente o sistema.
%                                             Deve ser um valor positivo. Se Ts <= 0 ou vazio, será usado Ts = 1
%                                             (frequência normalizada em rad/amostra).
%                                             (Padrão: 1)
%                                     - 'plot_results': Lógico (true/false ou 0/1). Se true, a rotina gera um gráfico
%                                                       da Norma H2 calculada para cada combinação alpha varrida.
%                                                       (Padrão: true)
%
% output:
%   output.max_h2_norm -> A norma H2 de pior caso encontrada em toda a varredura do simplex.
%   output.alpha_star  -> O vetor alpha (combinação convexa dos vértices) associado
%                         ao 'max_h2_norm'.
%
% Date: 25/09/2025
% Autor: glauco.gcps@gmail.com

%   Bloco de Validação de Entrada
if nargin < 4
    error('norma_h2_sistema_incerto_discrete:NotEnoughInputs', 'São necessários A_cell, B_cell, C_cell, D_cell.');
end

input_cells = {A_cell, B_cell, C_cell, D_cell};
input_names = {'A_cell', 'B_cell', 'C_cell', 'D_cell'};
num_vertices = numel(A_cell); % N, o número de vértices do politopo

% Valida que todas as entradas são cell arrays e têm o mesmo número de elementos
for i = 1:length(input_cells)
    if ~iscell(input_cells{i})
        error('norma_h2_sistema_incerto_discrete:InvalidInputType', 'O argumento %s deve ser um cell array.', input_names{i});
    end
    if numel(input_cells{i}) ~= num_vertices
        error('norma_h2_sistema_incerto_discrete:DimensionMismatch', 'Todos os cell arrays (A, B, C, D) devem ter o mesmo número de vértices (elementos).');
    end
end

% Valida as dimensões e tipos das matrizes dentro dos cell arrays
first_A_dims = size(A_cell{1});
first_B_dims = size(B_cell{1});
first_C_dims = size(C_cell{1});
first_D_dims = size(D_cell{1});

% Verifica a compatibilidade para o sistema de espaço de estados (A,B,C,D)
% A: n x n, B: n x m, C: p x n, D: p x m
if first_A_dims(1) ~= first_A_dims(2)
    error('norma_h2_sistema_incerto_discrete:InvalidADim', 'Matrizes A devem ser quadradas (n x n).');
end
if first_B_dims(1) ~= first_A_dims(1)
    error('norma_h2_sistema_incerto_discrete:IncompatibleDims', 'Dimensões de B (linhas) incompatíveis com A (linhas).');
end
if first_C_dims(2) ~= first_A_dims(1)
    error('norma_h2_sistema_incerto_discrete:IncompatibleDims', 'Dimensões de C (colunas) incompatíveis com A (colunas).');
end
if first_D_dims(1) ~= first_C_dims(1) || first_D_dims(2) ~= first_B_dims(2)
    error('norma_h2_sistema_incerto_discrete:IncompatibleDims', 'Dimensões de D incompatíveis com C e B.');
end

% Para a norma H2 discreta, D pode ser não-zero.
for k = 1:num_vertices
    current_A = A_cell{k};
    current_B = B_cell{k};
    current_C = C_cell{k};
    current_D = D_cell{k};

    if ~ismatrix(current_A) || ~isnumeric(current_A) || ~ismatrix(current_B) || ~isnumeric(current_B) || ...
       ~ismatrix(current_C) || ~isnumeric(current_C) || ~ismatrix(current_D) || ~isnumeric(current_D)
        error('norma_h2_sistema_incerto_discrete:InvalidMatrixElement', 'Todos os elementos dos cell arrays devem ser matrizes numéricas.');
    end
    
    if ~isequal(size(current_A), first_A_dims) || ~isequal(size(current_B), first_B_dims) || ...
       ~isequal(size(current_C), first_C_dims) || ~isequal(size(current_D), first_D_dims)
        error('norma_h2_sistema_incerto_discrete:DimensionMismatch', 'Todas as matrizes dentro de cada cell array devem ter as mesmas dimensões.');
    end
end
%   Fim do Bloco de Validação de Entrada  

%   Tratamento das Opções (varargin)  
options = [];
if nargin > 4 % A_cell, B_cell, C_cell, D_cell são 4 argumentos obrigatórios
    if nargin == 5 && isstruct(varargin{1}) % A, B, C, D, options_struct
        options = varargin{1};
    else % A, B, C, D, 'option_name', option_value, ...
        try
            options = struct(varargin{:});
        catch ME
            error('norma_h2_sistema_incerto_discrete:InvalidOptions', 'Argumentos opcionais inválidos. Use uma estrutura ou pares nome-valor. Erro: %s', ME.message);
        end
    end
end

%   Definição de Valores Padrão para as Opções  
if ~isfield(options, 'simplex_step')
    options.simplex_step = 0.001;
end
if ~isfield(options, 'Ts')
    options.Ts = 1; % Padrão: 1 (frequência normalizada)
end
if ~isfield(options, 'plot_results')
    options.plot_results = true; % Padrão: plotar os resultados
end

%   Validação das Opções  
if ~isnumeric(options.simplex_step) || ~isscalar(options.simplex_step) || options.simplex_step <= 0 || options.simplex_step > 1
    error('norma_h2_sistema_incerto_discrete:InvalidSimplexStep', 'simplex_step deve ser um escalar numérico positivo menor ou igual a 1.');
end
if ~isnumeric(options.Ts) || ~isscalar(options.Ts) || options.Ts <= 0
    warning('norma_h2_sistema_incerto_discrete:InvalidTs', 'Ts deve ser um escalar numérico positivo. Usando Ts = 1 (frequência normalizada).');
    options.Ts = 1; % Força para 1 para garantir natureza discreta
end
if ~islogical(options.plot_results) && ~ismember(options.plot_results, [0, 1])
    error('norma_h2_sistema_incerto_discrete:InvalidPlotResults', 'plot_results deve ser um valor lógico (true/false ou 0/1).');
end

%   Uso das Opções  
passo_simplex = options.simplex_step;
Ts = options.Ts;
plot_results = options.plot_results;

%   Geração dos pontos do Simplex  
% Assumimos que 'particao_simplex' está corretamente implementada e disponível.
alphas_simplex_points = particao_simplex(num_vertices, passo_simplex);

%   Inicialização para o Plot  
if plot_results
    figure;
    hold on;
    grid on;
    xlabel('Índice da Combinação Alpha'); % Eixo X para identificar a ordem da combinação alpha
    ylabel('Norma H2');
    title('Norma H2 para Diferentes Combinações de Alpha (Sistema Discreto)');
end

%   Variáveis para armazenar o máximo valor da norma H2 e o alpha associado  
max_h2_norm = -inf;
alpha_star = [];
h2_norms_per_alpha = zeros(1, length(alphas_simplex_points)); % Para armazenar todas as normas H2 para plotagem

%   Loop sobre os pontos do Simplex  
for k = 1:length(alphas_simplex_points)
    current_alpha_vector = alphas_simplex_points{k}; % Vetor coluna [alpha_1; ...; alpha_N]
    
    %   Construção das Matrizes Incertas (A,B,C,D)(alpha)  
    A_alpha = zeros(first_A_dims);
    B_alpha = zeros(first_B_dims);
    C_alpha = zeros(first_C_dims);
    D_alpha = zeros(first_D_dims);
    
    for i = 1:num_vertices
        A_alpha = A_alpha + current_alpha_vector(i) * A_cell{i};
        B_alpha = B_alpha + current_alpha_vector(i) * B_cell{i};
        C_alpha = C_alpha + current_alpha_vector(i) * C_cell{i};
        D_alpha = D_alpha + current_alpha_vector(i) * D_cell{i};
    end
    
    %   Criação do Objeto de Sistema de Espaço de Estados Discreto  
    % Passamos Ts para garantir que o sistema seja discreto.
    sys_alpha = ss(A_alpha, B_alpha, C_alpha, D_alpha, Ts);
    
    %   CÁLCULO DA NORMA H2  
    % A função norm(sys, 2) calcula a norma H2.
    % Para sistemas discretos, D pode ser não-zero.
    % A função norm() retornará 'inf' se o sistema não for estável.
    current_h2_norm = norm(sys_alpha, 2); 
    h2_norms_per_alpha(k) = current_h2_norm; % Armazena para plotagem
    
    %   Atualização do Máximo Valor da Norma H2  
    if current_h2_norm > max_h2_norm
        max_h2_norm = current_h2_norm;
        alpha_star = current_alpha_vector;
    end
end

%   Plot dos Resultados (se ativado)  
if plot_results
    plot(1:length(alphas_simplex_points), h2_norms_per_alpha, 'o', 'MarkerSize', 5, 'DisplayName', 'Norma H2 por Alpha');
    hold off;
end

%   Saída da Rotina  
output.max_h2_norm = max_h2_norm;
output.alpha_star = alpha_star;

fprintf('Processamento concluído para sistema incerto discreto.\n');
fprintf('Norma H2 de pior caso: %.4f\n', max_h2_norm);
fprintf('Alpha* associado: [%s]\n', num2str(alpha_star'));

end