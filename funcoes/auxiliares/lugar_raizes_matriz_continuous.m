function output = lugar_raizes_matriz_continuous(A, varargin)
% function output = lugar_raizes_matriz_continuous(A, varargin)
%
% Avalia o lugar das raízes de um sistema linear contínuo no tempo com incerteza politópica.
% A rotina varre o simplex unitário (combinações convexas dos vértices) para encontrar
% o vetor alpha que resulta no autovalor com a maior parte real.
%
% input:
%   A        -> cell array de matrizes de estado Ai (vértices do politopo).
%               Cada Ai deve ser uma matriz numérica quadrada e não vazia.
%               Todas as matrizes em A devem ter as mesmas dimensões (n x n, onde n é a ordem do sistema).
%   varargin -> (opcional) Estrutura de opções ou pares nome-valor para configurar a rotina.
%               Campos suportados:
%               - 'simplex_step': Scalar numérico. Define o passo para a geração de pontos no simplex unitário.
%                                 Um valor menor resulta em uma varredura mais fina e precisa do politopo,
%                                 mas aumenta significativamente o tempo de computação. Deve ser um valor
%                                 positivo e menor ou igual a 1.
%                                 (Padrão: 0.001)
%               - 'plot_results': Lógico (true/false ou 0/1). Se true, a rotina gera um gráfico do lugar das raízes
%                                 no plano complexo, incluindo os eixos real e imaginário.
%                                 (Padrão: true)
%
% output:
%   output.max_real_part -> O maior valor da parte real dos autovalores encontrados em toda a varredura do simplex.
%                           Este valor indica a estabilidade de pior caso dentro do politopo amostrado.
%   output.alpha_star    -> O vetor alpha (combinação convexa dos vértices) associado ao 'max_real_part'.
%                           Este vetor identifica a configuração específica do sistema que resulta na pior estabilidade.
%   output.LAMBDA        -> Um vetor coluna contendo todos os autovalores calculados durante a varredura do simplex.
%                           Pode ser usado para plotagens personalizadas ou análises adicionais.
%
% Date: 25/09/2025
% Autor: glauco.gcps@gmail.com

%   Bloco de Validação de Entrada 
if nargin < 1
    error('lugar_raizes_matriz_continuous:NotEnoughInputs', 'É necessário pelo menos um argumento de entrada (o cell array A).');
end

if ~iscell(A)
    error('lugar_raizes_matriz_continuous:InvalidInputType', 'O argumento de entrada A deve ser um cell array.');
end

if isempty(A)
    error('lugar_raizes_matriz_continuous:EmptyCellArray', 'O cell array A não pode estar vazio. Ele deve conter pelo menos uma matriz.');
end

num_matrices = numel(A); % N, o número de vértices do politopo
first_matrix_dims = [];

for i = 1:num_matrices
    current_matrix = A{i};
    
    if ~ismatrix(current_matrix) || isempty(current_matrix) || ~isnumeric(current_matrix)
        error('lugar_raizes_matriz_continuous:InvalidMatrixElement', ...
              'O elemento %d do cell array A não é uma matriz numérica válida e não vazia.', i);
    end
    
    if i == 1
        first_matrix_dims = size(current_matrix);
    else
        if ~isequal(size(current_matrix), first_matrix_dims)
            error('lugar_raizes_matriz_continuous:DimensionMismatch', ...
                  'As matrizes no cell array A devem ter as mesmas dimensões. Matriz %d tem dimensões diferentes.', i);
        end
    end
end
%   Fim do Bloco de Validação de Entrada  

%   Tratamento das Opções (varargin)  
options = [];
if nargin > 1
    if nargin == 2 % A, options_struct
        options = varargin{1};
        if ~isstruct(options)
            error('lugar_raizes_matriz_continuous:InvalidOptions', 'O segundo argumento deve ser uma estrutura de opções ou pares nome-valor.');
        end
    else % A, 'option_name', option_value, ...
        try
            options = struct(varargin{:});
        catch ME
            error('lugar_raizes_matriz_continuous:InvalidOptions', 'Argumentos opcionais inválidos. Use uma estrutura ou pares nome-valor. Erro: %s', ME.message);
        end
    end
end

%   Definição de Valores Padrão para as Opções  
if ~isfield(options, 'simplex_step')
    options.simplex_step = 0.001;
end
if ~isfield(options, 'plot_results')
    options.plot_results = true; % Padrão: plotar os resultados
end

%   Validação das Opções  
if ~isnumeric(options.simplex_step) || ~isscalar(options.simplex_step) || options.simplex_step <= 0 || options.simplex_step > 1
    error('lugar_raizes_matriz_continuous:InvalidSimplexStep', 'simplex_step deve ser um escalar numérico positivo menor ou igual a 1.');
end
if ~islogical(options.plot_results) && ~ismember(options.plot_results, [0, 1])
    error('lugar_raizes_matriz_continuous:InvalidPlotResults', 'plot_results deve ser um valor lógico (true/false ou 0/1).');
end

%   Uso das Opções  
passo_simplex = options.simplex_step;
plot_results = options.plot_results;

%   Geração dos pontos do Simplex  
% num_dimensoes para particao_simplex é o número de vértices (num_matrices)
% A função particao_simplex deve ser corrigida para retornar um cell array de vetores COLUNA
% Assumindo que particao_simplex é uma função auxiliar fornecida e corrigida para este formato.
alphas_simplex_points = particao_simplex(num_matrices, passo_simplex);

%   Variáveis para armazenar o máximo valor da parte real e o alpha associado  
max_real_part = -inf;
alpha_star = [];
LAMBDA = [];

%   Loop sobre os pontos do Simplex  
% alphas_simplex_points é um cell array de vetores coluna
for k = 1:length(alphas_simplex_points)
    current_alpha_vector = alphas_simplex_points{k}; % Este é o vetor [alpha_1; alpha_2; ...; alpha_N]
    
    %   Construção da Matriz Incerta A(alpha)  
    % A(alpha) = sum(alpha_i * Ai)
    A_incerta = zeros(first_matrix_dims); % Inicializa com zeros na dimensão correta
    
    for i = 1:num_matrices
        A_incerta = A_incerta + current_alpha_vector(i) * A{i};
    end
    
    %   Cálculo dos Autovalores  
    lambda_eigenvalues = eig(A_incerta);
    LAMBDA = [LAMBDA; lambda_eigenvalues]; 
    
    %   Atualização do Máximo Valor da Parte Real  
    current_max_real_part_for_alpha = max(real(lambda_eigenvalues));
        
    if current_max_real_part_for_alpha > max_real_part
        max_real_part = current_max_real_part_for_alpha;
        alpha_star = current_alpha_vector;
    end
end

%   Saída da Rotina  
output.max_real_part = max_real_part;
output.alpha_star = alpha_star;
output.LAMBDA = LAMBDA;

%   Plot dos Resultados (se ativado)  
if plot_results
    figure;
    hold on; % Permite plotar múltiplos pontos no mesmo gráfico
    grid on; % Adiciona grade para melhor visualização
    xlabel('Parte Real');
    ylabel('Parte Imaginária');
    title('Lugar das Raízes da Matriz Politópica A(\alpha)');
    
    % Plota todos os autovalores encontrados
    plot(real(LAMBDA), imag(LAMBDA), '.');
    
    % Adiciona os eixos real e imaginário para referência
    % Ajusta os limites para garantir que os eixos sejam visíveis
    min_real = min(real(LAMBDA));
    max_real = max(real(LAMBDA));
    min_imag = min(imag(LAMBDA));
    max_imag = max(imag(LAMBDA));
    
    plot([min_real - abs(min_real)*0.1 - 0.1, max_real + abs(max_real)*0.1 + 0.1], [0, 0], 'k--'); % Eixo real
    plot([0, 0], [min_imag - abs(min_imag)*0.1 - 0.1, max_imag + abs(max_imag)*0.1 + 0.1], 'k--'); % Eixo imaginário
    
    hold off;
end

fprintf('Processamento concluído. Máximo valor da parte real dos autovalores: %.4f\n', max_real_part);
fprintf('Alpha* associado: [%s]\n', num2str(alpha_star'));

end
