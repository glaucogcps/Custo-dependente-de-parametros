function output = hinf_lmi_mjls_d(A, E, Cz, Ez, Gamma, mu, varargin)
% function output = hinf_lmi_mjls_d(A,E,Cz,Ez,Gamma,varargin)
%
% Evaluate the H-infinity norm of a discrete-time MJLS. 
% Based on Lemma 1 for precisely known transition probabilities.
%
% input:  (A,E,Cz,Ez,Gamma) -> state-space matrices, probability matrix
%         Initial distribution 'mu'
%
% varargin (Opções):
%   .solver    : Solver LMI a ser utilizado (default: 'sedumi')
%   .hinf      : Limite superior desejado para a norma H-inf. Se for 0, a 
%                função minimiza a norma H-inf (default: 0)
%   .tolerance : Tolerância numérica para a verificação de viabilidade 
%                das LMIs (default: 1e-7)
%        
% output: output.hinf     -> H-infinity norm (0 if unfeasible)
%         output.cpusec   -> cpu time to solve the LMIs (seconds)
%         output.cpusec_m -> cpu time to mount the LMIs (seconds)
%         output.P        -> Lyapunov matrices
%         output.V        -> number of decision variables
%         output.L        -> number of LMIs rows
%
% Author: glauco.gcps@gmail.com
% Date: 03/04/2026

%   Bloco de Validação de Entrada
if nargin < 6
    error('hinf_mjls:NotEnoughInputs', 'São necessários A, E, Cz, Ez, Gamma e mu.');
end

input_cells = {A, E, Cz, Ez};
input_names = {'A', 'E', 'Cz', 'Ez'};
num_modos = numel(A); % N (Número de modos Markovianos / vértices)

% 1. Valida tipos e consistência do número de modos
for i = 1:length(input_cells)
    if ~iscell(input_cells{i})
        error('hinf_mjls:InvalidInputType', 'O argumento %s deve ser um cell array.', input_names{i});
    end
    if numel(input_cells{i}) ~= num_modos
        error('hinf_mjls:DimensionMismatch', 'Todos os cell arrays devem ter o mesmo número de modos (%d).', num_modos);
    end
end

% 2. Extrai dimensões baseando-se no primeiro modo
nx = size(A{1}, 1); % Número de estados
nw = size(E{1}, 2); % Número de entradas de distúrbio
nz = size(Cz{1}, 1); % Número de saídas de desempenho

dims_A  = size(A{1});
dims_E  = size(E{1});
dims_Cz = size(Cz{1});
dims_Ez = size(Ez{1});

% 3. Validações de coerência das dimensões (Álgebra Linear)
if dims_A(1) ~= dims_A(2)
    error('hinf_mjls:InvalidADim', 'Matrizes A devem ser quadradas (nx x nx).');
end
if dims_E(1) ~= nx
    error('hinf_mjls:IncompatibleDims', 'O número de linhas de E deve ser igual ao de A (nx = %d).', nx);
end
if dims_Cz(2) ~= nx
    error('hinf_mjls:IncompatibleDims', 'O número de colunas de Cz deve ser igual ao de A (nx = %d).', nx);
end
if dims_Ez(1) ~= nz || dims_Ez(2) ~= nw
    error('hinf_mjls:IncompatibleDims', 'Dimensões de Ez devem ser (nz x nw) -> (%d x %d).', nz, nw);
end

% 4. Valida conteúdo numérico e constância das dimensões em TODOS os modos
for k = 1:num_modos
    if ~isnumeric(A{k}) || ~isnumeric(E{k}) || ~isnumeric(Cz{k}) || ~isnumeric(Ez{k})
        error('hinf_mjls:NonNumeric', 'O conteúdo de todas as células no modo %d deve ser numérico.', k);
    end
    
    if ~isequal(size(A{k}), dims_A) || ~isequal(size(E{k}), dims_E) || ...
       ~isequal(size(Cz{k}), dims_Cz) || ~isequal(size(Ez{k}), dims_Ez)
        error('hinf_mjls:VaryingDims', 'As dimensões das matrizes não podem variar entre os modos (verifique o modo %d).', k);
    end
end

% 5. Validação das Propriedades Markovianas (Gamma e mu)
if size(Gamma, 1) ~= num_modos || size(Gamma, 2) ~= num_modos
    error('hinf_mjls:InvalidGamma', 'A matriz Gamma deve ser quadrada com dimensão igual ao número de modos (%d x %d).', num_modos, num_modos);
end
% Tolerância para soma das probabilidades devido a arredondamentos numéricos
if any(Gamma(:) < 0) || any(abs(sum(Gamma, 2) - 1) > 1e-6)
    error('hinf_mjls:InvalidGammaProb', 'Gamma deve ser uma matriz de transição estocástica válida (elementos >= 0 e soma das linhas = 1).');
end

if length(mu) ~= num_modos
    error('hinf_mjls:InvalidMu', 'O vetor mu deve ter comprimento igual ao número de modos (%d).', num_modos);
end
if any(mu < 0) || abs(sum(mu) - 1) > 1e-6
    error('hinf_mjls:InvalidMuProb', 'mu deve ser uma distribuição de probabilidade válida (elementos >= 0 e soma = 1).');
end

%   Tratamento das Opções (varargin)  
options = [];
if nargin > 6
    if nargin == 7 && isstruct(varargin{1})
        options = varargin{1};
    else
        try
            options = struct(varargin{:});
        catch
            error('hinf_mjls:InvalidOptions', 'Opções inválidas. Use struct ou pares nome-valor.');
        end
    end
end

if ~isfield(options,'solver'), options.solver = 'sedumi'; end
if ~isfield(options,'hinf'), options.hinf = 0; end
if ~isfield(options,'tolerance'), options.tolerance = 1e-7; end

%   Início da Lógica da Função (LMIs)
sigma = num_modos; 
output.cpusec_m = clock;

%new LMI system
LMIs = [];

%create the variables
if options.hinf == 0
    rho = sdpvar(1);
    obj = rho;
else
    rho = options.hinf * options.hinf;
    options.tolerance = 0;
    obj = [];
end

Xcal = [];
for i=1:sigma
    P{i} = sdpvar(nx,nx,'symmetric');
    Xcal = [Xcal; P{i}]; % Empilha as variáveis P
end

for i=1:sigma
    % Cria a matriz bloco Lambda usando Kronecker
    Lambda{i} = kron(Gamma(i,:), eye(nx));
    Pp{i} = Lambda{i} * Xcal;
end

for i=1:sigma
    Z11 = A{i};
    Z12 = E{i};
    Z21 = Cz{i};
    Z22 = Ez{i};
    Z = [Z11, Z12;
         Z21, Z22];
    
    X = Z';
    
    Y11 = Pp{i};
    Y12 = zeros(nx, nz);
    Y22 = eye(nz);
    Y = [Y11, Y12;
         Y12', Y22];
    
    W11 = P{i};
    W12 = zeros(nx, nw);
    W21 = zeros(nw, nx);
    W22 = rho*eye(nw);
    W = [W11, W12;
         W21, W22];
         
    LMIs = [LMIs, X*Y*Z - W <= 0];
end

output.cpusec_m = etime(clock,output.cpusec_m);

% determine the number of LMI rows
output.L = 0;
for i=1:size(LMIs,1)
    output.L = output.L + size(sdpvar(LMIs(i)),1);
end

% determine the number of scalar variables
output.V = size(getvariables(LMIs),2);

sol = solvesdp(LMIs, obj, sdpsettings('verbose',0,'solver',options.solver));
output.cpusec = sol.solvertime;

p = min(checkset(LMIs));
output.delta = p;
output.feas = 0;

%capturing the solutions (if ones exist)
if p > -options.tolerance
    output.feas = 1;
    output.hinf = sqrt(double(rho));   
    for i=1:sigma
        output.P{i} = double(P{i});
    end
end