function output = hinf_sistema_incerto_discrete(A_cell, B_cell, C_cell, D_cell, varargin)
% function output = valor_singular_sistema_incerto_discrete(A_cell, B_cell, C_cell, D_cell, varargin)
%
% Calcula o valor singular máximo (norma H_infinito) de um sistema incerto discreto no tempo
% (A,B,C,D)(alpha) e plota os diagramas de valor singular (magnitude).
%
% Esta função varre o simplex unitário para combinações de alpha, constrói
% o sistema (A,B,C,D)(alpha) e calcula o valor singular máximo (norma H_infinito)
% de pior caso, além de plotar os diagramas de valor singular (magnitude) para cada combinação.
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
%                                             Essencial para a função 'sigma' interpretar corretamente as frequências.
%                                             Deve ser um valor positivo. Se Ts <= 0 ou vazio, será usado Ts = 1
%                                             (frequência normalizada em rad/amostra).
%                                             (Padrão: 1)
%                                     - 'plot_legend': Lógico (true/false ou 0/1). Se true, a legenda do gráfico
%                                                      será exibida, identificando cada curva de valor singular.
%                                                      (Padrão: true)
%                                     - 'freq_range': Vetor de 2 elementos [Wmin, Wmax]. Define o intervalo de
%                                                     frequências (em rad/amostra) para o cálculo dos valores singulares.
%                                                     Para sistemas discretos, o intervalo válido é [0, pi/Ts] ou [0, pi]
%                                                     se Ts=1 (frequência normalizada).
%                                                     (Padrão: [0, pi])
%                                     - 'plot_graphic': Lógico (true/false ou 0/1). Se true, o gráfico de valor
%                                                       singular será plotado. Se false, a plotagem será suprimida.
%                                                       (Padrão: true)
%
% output:
%   output.max_h_inf_norm -> O valor singular máximo (norma H_infinito) de pior caso
%                            encontrado em toda a varredura do simplex.
%   output.alpha_star     -> O vetor alpha (combinação convexa dos vértices) associado
%                            ao 'max_h_inf_norm'.
%
% Date: 25/09/2025
% Autor: glauco.gcps@gmail.com

% Bloco de Validação de Entrada
if nargin < 4
    error('valor_singular_sistema_incerto_discrete:NotEnoughInputs', 'São necessários A_cell, B_cell, C_cell, D_cell.');
end

input_cells = {A_cell, B_cell, C_cell, D_cell};
input_names = {'A_cell', 'B_cell', 'C_cell', 'D_cell'};
num_vertices = numel(A_cell); % N, o número de vértices do politopo

% Valida que todas as entradas são cell arrays e têm o mesmo número de elementos
for i = 1:length(input_cells)
    if ~iscell(input_cells{i})
        error('valor_singular_sistema_incerto_discrete:InvalidInputType', 'O argumento %s deve ser um cell array.', input_names{i});
    end
    if numel(input_cells{i}) ~= num_vertices
        error('valor_singular_sistema_incerto_discrete:DimensionMismatch', 'Todos os cell arrays (A, B, C, D) devem ter o mesmo número de vértices (elementos).');
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
    error('valor_singular_sistema_incerto_discrete:InvalidADim', 'Matrizes A devem ser quadradas (n x n).');
end
if first_B_dims(1) ~= first_A_dims(1)
    error('valor_singular_sistema_incerto_discrete:IncompatibleDims', 'Dimensões de B (linhas) incompatíveis com A (linhas).');
end
if first_C_dims(2) ~= first_A_dims(1)
    error('valor_singular_sistema_incerto_discrete:IncompatibleDims', 'Dimensões de C (colunas) incompatíveis com A (colunas).');
end
if first_D_dims(1) ~= first_C_dims(1) || first_D_dims(2) ~= first_B_dims(2)
    error('valor_singular_sistema_incerto_discrete:IncompatibleDims', 'Dimensões de D (linhas/colunas) incompatíveis com C (linhas) e B (colunas).');
end

for k = 1:num_vertices
    current_A = A_cell{k};
    current_B = B_cell{k};
    current_C = C_cell{k};
    current_D = D_cell{k};

    if ~ismatrix(current_A) || ~isnumeric(current_A) || ~ismatrix(current_B) || ~isnumeric(current_B) || ...
       ~ismatrix(current_C) || ~isnumeric(current_C) || ~ismatrix(current_D) || ~isnumeric(current_D)
        error('valor_singular_sistema_incerto_discrete:InvalidMatrixElement', 'Todos os elementos dos cell arrays devem ser matrizes numéricas.');
    end
    
    if ~isequal(size(current_A), first_A_dims) || ~isequal(size(current_B), first_B_dims) || ...
       ~isequal(size(current_C), first_C_dims) || ~isequal(size(current_D), first_D_dims)
        error('valor_singular_sistema_incerto_discrete:DimensionMismatch', 'Todas as matrizes dentro de cada cell array devem ter as mesmas dimensões.');
    end
end
% Fim do Bloco de Validação de Entrada

% Tratamento das Opções (varargin)
options = [];
if nargin > 4 % A_cell, B_cell, C_cell, D_cell são 4 argumentos obrigatórios
    if nargin == 5 && isstruct(varargin{1}) % A, B, C, D, options_struct
        options = varargin{1};
    else % A, B, C, D, 'option_name', option_value, ...
        try
            options = struct(varargin{:});
        catch ME
            error('valor_singular_sistema_incerto_discrete:InvalidOptions', 'Argumentos opcionais inválidos. Use uma estrutura ou pares nome-valor. Erro: %s', ME.message);
        end
    end
end

% Definição de Valores Padrão para as Opções
if ~isfield(options, 'simplex_step')
    options.simplex_step = 0.001;
end
if ~isfield(options, 'Ts')
    options.Ts = 1; % Padrão: 1 (frequência normalizada)
end
if ~isfield(options, 'plot_legend')
    options.plot_legend = true; % Padrão: exibir legenda
end
if ~isfield(options, 'freq_range')
    options.freq_range = [0, pi]; % Padrão: de 0 a pi rad/amostra (para Ts=1)
end
if ~isfield(options, 'plot_graphic')
    options.plot_graphic = true; % Padrão: plotar gráfico por padrão
end
if isfield(options, 'hinf')
    hinf_incerto = true; 
    hinf = options.hinf;
else
    hinf_incerto = false;
end

% Validação das Opções
if ~isnumeric(options.simplex_step) || ~isscalar(options.simplex_step) || options.simplex_step <= 0 || options.simplex_step > 1
    error('valor_singular_sistema_incerto_discrete:InvalidSimplexStep', 'simplex_step deve ser um escalar numérico positivo menor ou igual a 1.');
end
if ~isnumeric(options.Ts) || ~isscalar(options.Ts) || options.Ts <= 0
    warning('valor_singular_sistema_incerto_discrete:InvalidTs', 'Ts deve ser um escalar numérico positivo. Usando Ts = 1 (frequência normalizada).');
    options.Ts = 1; % Força para 1 para garantir natureza discreta
end
if ~islogical(options.plot_legend) && ~ismember(options.plot_legend, [0, 1])
    error('valor_singular_sistema_incerto_discrete:InvalidPlotLegend', 'plot_legend deve ser um valor lógico (true/false ou 0/1).');
end
if ~isnumeric(options.freq_range) || ~isvector(options.freq_range) || numel(options.freq_range) ~= 2 || options.freq_range(1) < 0 || options.freq_range(2) <= options.freq_range(1)
    error('valor_singular_sistema_incerto_discrete:InvalidFreqRange', 'freq_range deve ser um vetor [Wmin, Wmax] com Wmin >= 0 e Wmax > Wmin.');
end
% Validação específica para freq_range em sistemas discretos: max_freq <= pi/Ts
if options.freq_range(2) > pi/options.Ts + 1e-9 % Adiciona tolerância para ponto flutuante
    warning('valor_singular_sistema_incerto_discrete:FreqRangeExceedsNyquist', 'Wmax excede a frequência de Nyquist (pi/Ts). Limitando Wmax a pi/Ts.');
    options.freq_range(2) = pi/options.Ts;
end
if ~islogical(options.plot_graphic) && ~ismember(options.plot_graphic, [0, 1])
    error('valor_singular_sistema_incerto_discrete:InvalidPlotGraphic', 'plot_graphic deve ser um valor lógico (true/false ou 0/1).');
end


%   Uso das Opções  
passo_simplex = options.simplex_step;
Ts = options.Ts;
plot_legend = options.plot_legend;
Wmin = options.freq_range(1);
Wmax = options.freq_range(2);
plot_graphic = options.plot_graphic; % NOVO: flag para plotagem

alphas_simplex_points = particao_simplex(num_vertices, passo_simplex);

%   Inicialização para o Plot (CONDICIONAL)  
if plot_graphic
    figure;
    hold on;
    grid on;
    xlabel('Frequência (rad/amostra)'); % Frequência para sistemas discretos
    ylabel('Magnitude (dB)');
    title('Diagrama de Valor Singular (Magnitude) para Sistema Incerto Discreto');
end

%   Variáveis para armazenar o máximo valor singular e o alpha associado  
max_h_inf_norm = -inf;
alpha_star = [];

% alpha, norma de pior caso, custo garantido
tab = [];

%   Loop sobre os pontos do Simplex  
% Coleta os handles das linhas e os nomes para a legenda (CONDICIONAL)
line_handles = gobjects(0); % Inicializa um array vazio de handles de objetos gráficos
legend_entries = {}; % Inicializa um cell array vazio para as entradas da legenda

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
    % O período de amostragem Ts é passado para indicar que é um sistema discreto
    sys_alpha = ss(A_alpha, B_alpha, C_alpha, D_alpha, Ts);
    
    %   Cálculo dos Valores Singulares (Magnitude)  
    % Para sistemas discretos, sigma() sem range varre de 0 a pi/Ts.
    % Se 'freq_range' for especificado, sigma() usa esse range.
    [mag, w] = sigma(sys_alpha, {Wmin, Wmax}); 
    
    %   Plot do Diagrama de Valor Singular de Magnitude (CONDICIONAL)  
    if plot_graphic
        % Armazena o handle da linha para uso posterior na legenda
        h_line = semilogx(w, 20*log10(mag(1,:))); 
        line_handles(end+1) = h_line; % Adiciona o handle ao array
        legend_entries{end+1} = sprintf('\alpha=[%s]', num2str(current_alpha_vector', 2)); % Adiciona a entrada da legenda
    end
    
    %   Atualização do Máximo Valor Singular (Norma H_infinito)  
    current_h_inf_norm = max(mag(1,:)); 
    tab = [tab; current_alpha_vector', current_h_inf_norm];
    if hinf_incerto
       mu = sqrt(evalpar(hinf,{current_alpha_vector}));
       
       tab = [tab; current_alpha_vector', current_h_inf_norm, mu]; 
    end
    if current_h_inf_norm > max_h_inf_norm
        max_h_inf_norm = current_h_inf_norm;
        alpha_star = current_alpha_vector;
    end
end

%   Finalização do Plot (CONDICIONAL)  
if plot_graphic
    drawnow; % Força a atualização do gráfico
    hold off; % Libera o gráfico
    %   Ativação/Desativação da Legenda com base na flag  
    if plot_legend
        legend(line_handles, legend_entries, 'Location', 'best'); 
    else
        legend('off'); 
    end
end

%   Saída da Rotina  
output.tab = tab;
output.max_h_inf_norm = max_h_inf_norm;
output.alpha_star = alpha_star;

fprintf('Processamento concluído para sistema incerto discreto.\n');
fprintf('Máximo valor singular (Norma H_infinito) de pior caso: %.4f (%.2f dB)\n', max_h_inf_norm, 20*log10(max_h_inf_norm));
fprintf('Alpha* associado: [%s]\n', num2str(alpha_star'));

end