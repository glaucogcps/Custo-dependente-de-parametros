function output = valor_singular_sistema_incerto(A_cell, B_cell, C_cell, D_cell, type_analysis, varargin)
% Calcula o Valor singular de uma matriz politópica
% para sistemas contínuos ou discretos.
%
% Esta função atua como um dispatcher, chamando a função apropriada
% (valor_singular_sistema_incerto_continuous ou valor_singular_sistema_incerto_discrete)
% com base no 'type_analysis' fornecido.
%
% Sintaxe:
%   output = valor_singular_sistema_incerto(A_cell, B_cell, C_cell, D_cell, type_analysis, optional_simplex_step)
%
% Entradas:
%   A_cell        - Um cell array de matrizes A (A_cell{1}, A_cell{2}, ..., A_cell{N})
%                   representando os vértices do politopo para a matriz de estados.
%                   Todas as matrizes devem ter as mesmas dimensões.
%   B_cell        - Um cell array de matrizes B (B_cell{1}, B_cell{2}, ..., B_cell{N})
%                   representando os vértices do politopo para a matriz de entrada.
%   C_cell        - Um cell array de matrizes C (C_cell{1}, C_cell{2}, ..., C_cell{N})
%                   representando os vértices do politopo para a matriz de saída.
%   D_cell        - Um cell array de matrizes D (D_cell{1}, D_cell{2}, ..., D_cell{N})
%                   representando os vértices do politopo para a matriz de feedthrough.
%   type_analysis - Uma string que especifica o tipo de análise:
%                   - 'continuous' ou 'c' para sistemas contínuos.
%                   - 'discrete' ou 'd' para sistemas discretos.
%   varargin      -> (opcional) Estrutura de opções ou pares nome-valor para configurar a rotina.
%                    Campos suportados:
%                    - 'simplex_step': Scalar numérico. Define o passo para a geração de pontos
%                                      no simplex unitário. Um valor menor resulta em uma varredura
%                                      mais fina e precisa do politopo, mas aumenta o tempo de computação.
%                                      Deve ser um valor positivo e menor ou igual a 1.
%                                      (Padrão: 0.001)
%                     - 'plot_legend': Lógico (true/false ou 0/1). Se true, a legenda do gráfico
%                                      será exibida, identificando cada curva de valor singular.
%                                      (Padrão: true)
%                     - 'freq_range': Vetor de 2 elementos [Wmin, Wmax]. Define o intervalo de
%                                     frequências (em rad/s) para o cálculo dos valores singulares.
%                                     (Padrão: [0, 100])
%                     - 'Ts': Scalar numérico. Período de amostragem do sistema discreto.
%                             Essencial para a função 'sigma' interpretar corretamente as frequências.
%                             Deve ser um valor positivo. Se Ts <= 0 ou vazio, será usado Ts = 1
%                             (frequência normalizada em rad/amostra).
%                             (Padrão: 1)
%                     - 'plot_graphic': Lógico (true/false ou 0/1). Se true, o gráfico de valor
%                                      singular será plotado. Se false, a plotagem será suprimida.
%                                      (Padrão: true)
%
% Saída:
%   output        - Uma estrutura contendo os resultados da análise,
%                   dependendo da função específica chamada.
%                   - Para 'continuous': .max_h_inf_norm, .alpha_star
%                   - Para 'discrete':   .max_h_inf_norm, .alpha_star
% Date: 25/09/2025
% Autor: glauco.gcps@gmail.com

%   1. Validação de Entrada para type_analysis  
if nargin < 5 % A_cell, B_cell, C_cell, D_cell, type_analysis
    error('valor_singular_sistema_incerto:NotEnoughInputs', ...
          'É necessário especificar as cell arrays de matrizes (A, B, C, D) e o tipo de análise (e.g., ''continuous'' ou ''discrete'').');
end

if ~ischar(type_analysis)
    error('valor_singular_sistema_incerto:InvalidTypeArgument', ...
          'O quinto argumento (tipo de análise) deve ser uma string.');
end

% Converte para minúsculas para garantir a insensibilidade a maiúsculas/minúsculas
type_analysis_lower = lower(type_analysis);

%   2. Dispatching para a Função Apropriada  
if strcmp(type_analysis_lower, 'continuous') || strcmp(type_analysis_lower, 'c')
    fprintf('Chamando valor_singular_sistema_incerto_continuous para análise de sistema contínuo.\n');
    output = valor_singular_sistema_incerto_continuous(A_cell, B_cell, C_cell, D_cell, varargin{:}); % Usa {:} para desempacotar varargin
elseif strcmp(type_analysis_lower, 'discrete') || strcmp(type_analysis_lower, 'd')
    fprintf('Chamando valor_singular_sistema_incerto_discrete para análise de sistema discreto.\n');
    output = valor_singular_sistema_incerto_discrete(A_cell, B_cell, C_cell, D_cell, varargin{:}); % Usa {:} para desempacotar varargin
else
    error('valor_singular_sistema_incerto:UnknownAnalysisType', ...
          'Tipo de análise desconhecido. Use ''continuous''/''c'' ou ''discrete''/''d''.');
end

end