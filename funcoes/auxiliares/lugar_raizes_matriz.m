function output = lugar_raizes_matriz(A, type_analysis, varargin)
% LUGAR_RAIZES_MATRIZ Calcula o lugar das raízes de uma matriz politópica
% para sistemas contínuos ou discretos.
%
% Esta função atua como um dispatcher, chamando a função apropriada
% (lugar_raizes_matriz_continuous ou lugar_raizes_matriz_discrete)
% com base no 'type_analysis' fornecido.
%
% Sintaxe:
%   output = lugar_raizes_matriz(A, type_analysis, optional_simplex_step)
%
% Entradas:
%   A             - Um cell array de matrizes (A{1}, A{2}, ..., A{N})
%                   representando os vértices do politopo. Todas as matrizes
%                   devem ter as mesmas dimensões.
%   type_analysis - Uma string que especifica o tipo de análise:
%                   - 'continuous' ou 'c' para sistemas contínuos.
%                   - 'discrete' ou 'd' para sistemas discretos.
%   varargin      - (opcional) Estrutura de opções ou pares nome-valor para configurar a rotina.
%               Campos suportados:
%               - 'simplex_step': Scalar numérico. Define o passo para a geração de pontos no simplex unitário.
%                                 Um valor menor resulta em uma varredura mais fina e precisa do politopo,
%                                 mas aumenta significativamente o tempo de computação. Deve ser um valor
%                                 positivo e menor ou igual a 1.
%                                 (Padrão: 0.001)
%               - 'plot_results': Lógico (true/false ou 0/1). Se true, a rotina gera um gráfico do lugar das raízes
%                                 no plano complexo, incluindo os eixos real e imaginário.
%                                 (Padrão: true).
%
% Saída:
%   output        - Uma estrutura contendo os resultados da análise,
%                   dependendo da função específica chamada.
%                   - Para 'continuous': .max_real_part, .alpha_star, .LAMBDA
%                   - Para 'discrete':   .max_abs_value, .alpha_star, .LAMBDA
% Data: 25/09/2025
% Autor: glauco.gcps@gmail.com

%   1. Validação de Entrada para type_analysis  
if nargin < 2
    error('lugar_raizes_matriz:NotEnoughInputs', ...
          'É necessário especificar o tipo de análise (e.g., ''continuous'' ou ''discrete'').');
end

if ~ischar(type_analysis)
    error('lugar_raizes_matriz:InvalidTypeArgument', ...
          'O segundo argumento (tipo de análise) deve ser uma string.');
end

% Converte para minúsculas para garantir a insensibilidade a maiúsculas/minúsculas
type_analysis_lower = lower(type_analysis);

%   2. Dispatching para a Função Apropriada  
if strcmp(type_analysis_lower, 'continuous') || strcmp(type_analysis_lower, 'c')
    fprintf('Chamando lugar_raizes_matriz_continuous para análise de sistema contínuo.\n');
    output = lugar_raizes_matriz_continuous(A, varargin{:}); % Usa {:} para desempacotar varargin
elseif strcmp(type_analysis_lower, 'discrete') || strcmp(type_analysis_lower, 'd')
    fprintf('Chamando lugar_raizes_matriz_discrete para análise de sistema discreto.\n');
    output = lugar_raizes_matriz_discrete(A, varargin{:}); % Usa {:} para desempacotar varargin
else
    error('lugar_raizes_matriz:UnknownAnalysisType', ...
          'Tipo de análise desconhecido. Use ''continuous''/''c'' ou ''discrete''/''d''.');
end

end
    