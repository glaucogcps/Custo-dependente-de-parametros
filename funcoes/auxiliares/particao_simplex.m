function pontos_simplex = particao_simplex(num_dimensoes, passo)
% particao_simplex Gera pontos de um simplex unitário com uma dada granularidade.
%
%   pontos_simplex = particao_simplex(num_dimensoes, passo)
%
%   Esta função gera um conjunto de pontos que representam uma partição
%   discreta de um simplex unitário (onde a soma das coordenadas é 1).
%   A granularidade da partição é definida pelo 'passo'.
%
%   Entradas:
%     num_dimensoes - Número de dimensões do simplex (inteiro positivo).
%     passo         - O tamanho do passo para a discretização das coordenadas
%                     (número positivo, e.g., 0.1, 0.05).
%
%   Saída:
%     pontos_simplex - Um cell array, onde cada célula contém um vetor coluna
%                      representando um ponto no simplex.
%
%   Otimização: A função verifica se um arquivo de cache (.mat) já existe
%   para a combinação específica de 'num_dimensoes' e 'passo'. Se existir,
%   carrega os dados do arquivo para evitar recálculos. Caso contrário,
%   realiza o cálculo completo e salva o resultado para uso futuro.
%
%   Esta versão utiliza um algoritmo ITERATIVO para otimizar o consumo de
%   memória e tempo em comparação com a versão recursiva.
% Date: 25/09/2025
% Autor: glauco.gcps@gmail.com

% Validação de entrada
if nargin < 2
    error('particao_simplex:NotEnoughInputs', 'São necessários num_dimensoes e passo.');
end
if ~isnumeric(num_dimensoes) || num_dimensoes < 1 || mod(num_dimensoes, 1) ~= 0
    error('particao_simplex:InvalidDims', 'num_dimensoes deve ser um inteiro positivo.');
end
if ~isnumeric(passo) || passo <= 0
    error('particao_simplex:InvalidStep', 'passo deve ser um número positivo.');
end

%   Geração do nome do arquivo de cache  
% Converte 'passo' para string, substituindo o ponto decimal por 'p' para o nome do arquivo
% Usamos '%.10f' para garantir precisão suficiente na representação do passo no nome do arquivo.
passo_str = strrep(num2str(passo, '%.10f'), '.', 'p'); 
filename = sprintf('alfas_simplex_dimensao_%d_passo_%s.mat', num_dimensoes, passo_str);

%   Verifica se o arquivo de cache existe  
if exist(filename, 'file') == 2
    fprintf('Arquivo de cache "%s" encontrado. Carregando dados...\n', filename);
    
    % Carrega o arquivo MAT para uma estrutura temporária
    temp_data = load(filename); 
    
    % Verifica se a variável esperada existe na estrutura
    if isfield(temp_data, 'pontos_simplex_saved')
        pontos_simplex = temp_data.pontos_simplex_saved; 
        fprintf('Dados carregados com sucesso do cache.\n');
        return; 
    else
        fprintf('AVISO: Variável "pontos_simplex_saved" não encontrada no cache. Recalculando...\n');
        % Se a variável não for encontrada, o código continua para recalcular
    end
else
    fprintf('Arquivo de cache "%s" não encontrado. Realizando cálculo (ITERATIVO)...\n', filename);
end

%   Início do Cálculo Iterativo (se o arquivo de cache não foi encontrado)  

soma_total = 1;
tolerance_float = 1e-9; % Tolerância para operações de ponto flutuante

% current_level_points armazena células no formato {vetor_prefixo, soma_restante}
% Começamos com um prefixo vazio e a soma total restante
current_level_points = {{[], soma_total}}; 

% Itera sobre cada dimensão para construir os pontos
for dim_idx = 1:num_dimensoes
    next_level_points = {}; % Lista para os pontos da próxima dimensão
    
    for k = 1:length(current_level_points)
        partial_vec = current_level_points{k}{1};
        remaining_sum = current_level_points{k}{2};
        
        % Se for a última dimensão, o valor é fixo (a soma restante)
        if dim_idx == num_dimensoes
            val = remaining_sum;
            val = round(val / passo) * passo; % Arredonda para múltiplo do passo
            
            % Trata valores negativos ou muito pequenos como zero
            if val < 0, val = 0; end
            if abs(val) < tolerance_float, val = 0; end
            
            new_point = [partial_vec; val];
            next_level_points{end+1} = {new_point, 0}; % A soma restante é 0
        else
            % Itera sobre os possíveis valores para a dimensão atual
            % Adiciona passo/2 para garantir que 'remaining_sum' seja incluído no loop
            for val = 0:passo:remaining_sum + passo/2 
                val = round(val / passo) * passo; % Arredonda para múltiplo do passo
                
                % Trata valores negativos ou muito pequenos como zero
                if val < 0, val = 0; end
                if abs(val) < tolerance_float, val = 0; end
                
                % Se o valor exceder a soma restante (devido a arredondamento), pula
                if val > remaining_sum + tolerance_float, continue; end
                
                new_partial_vec = [partial_vec; val];
                new_remaining_sum = remaining_sum - val;
                
                next_level_points{end+1} = {new_partial_vec, new_remaining_sum};
            end
        end
    end
    current_level_points = next_level_points; % Atualiza para a próxima iteração
end

%   Pós-processamento: Extrair pontos finais, aplicar correção de soma, filtrar e remover duplicatas  
final_points_matrix_rows = [];
tolerance_sum_check = 1e-9; % Tolerância para verificar a soma final

for k = 1:length(current_level_points)
    full_point_vec = current_level_points{k}{1}; % O vetor completo de N dimensões
    
    % Aplica a correção crucial para garantir que a soma seja exatamente 1
    % Isso é vital para a precisão e para evitar duplicatas por erro de ponto flutuante
    sum_of_first_n_minus_1_dims = sum(full_point_vec(1:end-1));
    full_point_vec(end) = soma_total - sum_of_first_n_minus_1_dims;
    
    % Arredonda novamente todos os elementos para garantir que são múltiplos de 'passo'
    full_point_vec = round(full_point_vec / passo) * passo;

    % Verifica se a soma do ponto (após o ajuste) está dentro da tolerância esperada
    if abs(sum(full_point_vec) - soma_total) < tolerance_sum_check
        final_points_matrix_rows = [final_points_matrix_rows; full_point_vec']; 
    end
end

% Remover duplicatas e garantir uma ordem consistente
if ~isempty(final_points_matrix_rows)
    unique_matrix_points_rows = unique(final_points_matrix_rows, 'rows');
    
    % Converte de volta para cell array de vetores COLUNA (formato de saída desejado)
    pontos_simplex = cell(1, size(unique_matrix_points_rows, 1)); 
    for k = 1:size(unique_matrix_points_rows, 1)
        pontos_simplex{k} = unique_matrix_points_rows(k, :)'; 
    end
else
    pontos_simplex = {}; 
end

%   Salva o resultado no arquivo de cache para uso futuro  
% Garante que apenas a variável 'pontos_simplex_saved' seja salva
pontos_simplex_saved = pontos_simplex; 
save(filename, 'pontos_simplex_saved'); 
fprintf('Resultado calculado e salvo em "%s".\n', filename);

end