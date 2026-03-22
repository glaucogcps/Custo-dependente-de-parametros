function pontos_simplex = particao_simplex_recursiva(num_dimensoes, passo)
% function pontos_simplex = particao_simplex_recursiva(num_dimensoes, passo)
% Algoritmo 2 desenvolvido
%
% Gera pontos de um simplex unitário com uma dada granularidade utilizando
% uma abordagem RECURSIVA.
%
% Esta função atua como um "wrapper" (casca). Ela chama a função auxiliar
% 'generate_simplex_points' para obter os pontos brutos e, em seguida,
% aplica um rigoroso pós-processamento para garantir que a soma das 
% coordenadas seja exatamente 1, eliminando erros de ponto flutuante e 
% removendo eventuais duplicatas.
%
% Entradas:
%   num_dimensoes - Número de dimensões do simplex (inteiro positivo).
%                   Representa o número de vértices do politopo (N).
%   passo         - O tamanho do passo para a discretização das coordenadas
%                   (número positivo, e.g., 0.1, 0.05).
%
% Saída:
%   pontos_simplex - Um cell array, onde cada célula contém um vetor coluna
%                    representando um ponto válido no simplex unitário.
%
% Date: 25/09/2025
% Autor: glauco.gcps@gmail.com

%   Validação de entrada  
if nargin < 2
    error('particao_simplex_recursiva:NotEnoughInputs', 'São necessários num_dimensoes e passo.');
end
if ~isnumeric(num_dimensoes) || num_dimensoes < 1 || mod(num_dimensoes, 1) ~= 0
    error('particao_simplex_recursiva:InvalidDims', 'num_dimensoes deve ser um inteiro positivo.');
end
if ~isnumeric(passo) || passo <= 0
    error('particao_simplex_recursiva:InvalidStep', 'passo deve ser um número positivo.');
end

%   1. Geração dos Pontos Brutos (Chamada Recursiva)  
soma_total = 1;
% Chama o motor recursivo iniciando com prefixo vazio
pontos_brutos = generate_simplex_points(soma_total, num_dimensoes, passo, []);

%   2. Pós-processamento e Correção Numérica  
% Inicializa a matriz para armazenar os pontos como linhas
M = zeros(length(pontos_brutos), num_dimensoes);

for i = 1:length(pontos_brutos)
    p = pontos_brutos{i}; % Vetor coluna bruto
    
    % Arredondamento inicial para alinhar à malha
    p = round(p / passo) * passo;
    
    % Correção CRUCIAL da última dimensão:
    % Garante que a soma de todas as coordenadas seja matematicamente igual a 1
    soma_n_menos_1 = sum(p(1:end-1));
    p(end) = soma_total - soma_n_menos_1;
    
    % Arredondamento final para garantir que a correção não tirou o ponto da malha
    p = round(p / passo) * passo;
    
    % Adiciona o vetor transposto como uma linha na matriz M
    M(i, :) = p';
end

%   3. Remoção de Duplicatas e Formatação de Saída  
if ~isempty(M)
    % Remove linhas duplicadas que podem ter surgido após a correção de arredondamento
    M_unico = unique(M, 'rows');
    
    % Converte a matriz única de volta para o formato de saída desejado (cell array de colunas)
    pontos_simplex = cell(1, size(M_unico, 1));
    for k = 1:size(M_unico, 1)
        pontos_simplex{k} = M_unico(k, :)';
    end
else
    pontos_simplex = {};
end

end