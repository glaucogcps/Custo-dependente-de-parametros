function all_points = generate_simplex_points_recursivo(current_sum, remaining_dims, step_size, current_prefix)

% Inicializa a lista de pontos a serem retornados
all_points = {}; 

% Caso Base: Se resta apenas 1 dimensão para preencher
if remaining_dims == 1
    % A última dimensão deve ser igual à soma restante
    % Arredondar para lidar com erros de ponto flutuante
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
    all_points{1} = [current_prefix; last_val]; % <--- MUDANÇA: Usando ';' para empilhar verticalmente
    return;
end

% Passo Recursivo: Itera sobre os possíveis valores para a próxima dimensão
% O valor pode ir de 0 até current_sum, em passos de step_size
% Usamos uma pequena tolerância para o loop for para evitar problemas de ponto flutuante
for val = 0:step_size:current_sum + step_size/2 % Adiciona step_size/2 para incluir o current_sum
    
    % Arredonda o valor para o múltiplo mais próximo de step_size
    val = round(val / step_size) * step_size;
    
    % Garante que o valor não seja negativo devido a erros de ponto flutuante
    if val < 0
        val = 0;
    end

    % Se o valor é muito próximo de zero, considere-o zero
    if abs(val) < 1e-9
        val = 0;
    end

    % Se o valor é maior que a soma restante, não é válido
    if val > current_sum + 1e-9 % Tolerância para current_sum
        continue; 
    end

    % Calcula a soma restante para as próximas dimensões
    next_sum = current_sum - val;
    
    % Constrói o novo prefixo para a próxima chamada recursiva como um VETOR COLUNA
    new_prefix = [current_prefix; val]; % <--- MUDANÇA: Usando ';' para empilhar verticalmente
    
    % Chama a função recursivamente para as dimensões restantes
    sub_points = generate_simplex_points(next_sum, remaining_dims - 1, step_size, new_prefix);
    
    % Adiciona os pontos gerados pela chamada recursiva à lista total
    all_points = [all_points, sub_points];
end

end