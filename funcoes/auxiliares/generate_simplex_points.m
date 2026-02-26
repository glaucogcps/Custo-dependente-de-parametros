function all_points = generate_simplex_points(current_sum, remaining_dims, step_size, current_prefix)
% function all_points = generate_simplex_points(current_sum, remaining_dims, step_size, current_prefix)
% Algoritmo 2
%
% Função auxiliar recursiva para gerar todos os pontos discretos em um simplex unitário.
% Os pontos são gerados de forma que a soma de suas coordenadas seja igual a 1,
% com cada coordenada sendo um múltiplo de 'step_size'.
%
% Esta função é tipicamente chamada por uma função wrapper (como 'particao_simplex')
% que inicializa 'current_sum' como 1 e 'current_prefix' como vazio.
%
% input:
%   current_sum    -> Scalar numérico. A soma que as 'remaining_dims' devem totalizar.
%                     Na chamada inicial para um simplex unitário, este valor é 1.
%   remaining_dims -> Scalar inteiro. O número de dimensões restantes a serem preenchidas
%                     para formar um ponto completo no simplex.
%   step_size      -> Scalar numérico. O incremento entre os valores possíveis para cada
%                     coordenada. Define a "malha" de discretização do simplex.
%   current_prefix -> Vetor coluna numérico. O prefixo parcial do ponto do simplex
%                     que está sendo construído. Contém os valores já determinados
%                     para as dimensões anteriores.
%
% output:
%   all_points     -> Cell array. Cada célula contém um vetor coluna numérico que
%                     representa um ponto completo no simplex unitário.
%                     Ex: {[alpha1; alpha2; ...; alphaN], [beta1; beta2; ...; betaN], ...}
%
% Date: 25/09/2025
% Autor: glauco.gcps@gmail.com

% Inicializa a lista de pontos a serem retornados
all_points = {}; 

% Caso Base: Se resta apenas 1 dimensão para preencher
if remaining_dims == 1
    % A última dimensão deve ser igual à soma restante
    % Arredondar para lidar com erros de ponto flutuante, garantindo que o valor
    % seja um múltiplo exato de step_size.
    last_val = round(current_sum / step_size) * step_size; 
    
    % Garante que o valor não seja negativo devido a erros de ponto flutuante
    if last_val < 0 
        last_val = 0; 
    end

    % Se o valor final é muito próximo de zero, considere-o zero
    if abs(last_val) < 1e-9 % Tolerância para zero
        last_val = 0;
    end

    % Se o valor final é muito próximo da soma restante, considere-o a soma
    if abs(last_val - current_sum) < 1e-9
        last_val = current_sum;
    end

    % Adiciona o ponto completo (prefixo + último valor) como um VETOR COLUNA
    all_points{1} = [current_prefix; last_val]; 
    return;
end

% Passo Recursivo: Itera sobre os possíveis valores para a próxima dimensão
% O valor pode ir de 0 até current_sum, em passos de step_size.
% Uma pequena tolerância (step_size/2) é adicionada ao limite superior do loop
% para garantir a inclusão do 'current_sum' devido a possíveis imprecisões de ponto flutuante.
for val = 0:step_size:current_sum + step_size/2 
    
    % Arredonda o valor para o múltiplo mais próximo de step_size para corrigir
    % eventuais erros de ponto flutuante acumulados.
    val = round(val / step_size) * step_size;
    
    % Garante que o valor não seja negativo devido a erros de ponto flutuante
    if val < 0
        val = 0;
    end

    % Se o valor é muito próximo de zero, considere-o zero
    if abs(val) < 1e-9
        val = 0;
    end

    % Se o valor é maior que a soma restante, não é válido para esta ramificação.
    % Uma pequena tolerância (1e-9) é usada para comparação.
    if val > current_sum + 1e-9 
        continue; 
    end

    % Calcula a soma restante para as próximas dimensões
    next_sum = current_sum - val;
    
    % Constrói o novo prefixo para a próxima chamada recursiva como um VETOR COLUNA
    new_prefix = [current_prefix; val]; 
    
    % Chama a função recursivamente para as dimensões restantes
    sub_points = generate_simplex_points(next_sum, remaining_dims - 1, step_size, new_prefix);
    
    % Adiciona os pontos gerados pela chamada recursiva à lista total
    all_points = [all_points, sub_points];
end

end
