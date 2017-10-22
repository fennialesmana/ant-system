clear;clc;
% INITIALIZE DATA
%nodes = csvread('city10.csv')';
%nodes = csvread('eil51.csv')';
nodes = csvread('berlin52.csv')';

% PARAMETER SETTING
MAX_ITERATION = 5000;
Q = 10;
alpha = 1;
beta = 2;
rho = 0.5;
%initial_tao = 0.001;

iteration = 1;
number_of_nodes = size(nodes, 2);
number_of_ants = number_of_nodes;
initial_tao = 1/(number_of_nodes*7542); % change it according to the dataset

tao = ones(number_of_nodes, number_of_nodes) * initial_tao;
distance = dist(nodes);
visibility = 1./distance;

shortest_path = zeros(1, number_of_nodes+1);
shortest_distance_each_iteration = zeros(MAX_ITERATION, 1);
shortest_distance_all_iteration = zeros(MAX_ITERATION, 1);

while iteration <= MAX_ITERATION
    s=1; % s: tabu list index -> s-th city will be filled
    tabu = zeros(number_of_ants, number_of_nodes+1); % save the route of each ant
    
    % random initial node of ant
    tabu(:, 1) = floor((rand(number_of_ants, 1)*number_of_nodes)+1);

    % fill tabu
    while s <= number_of_nodes
        s = s + 1; % s = tabu index to be filled
        for k = 1:number_of_ants % looping according to total ants | canculate the probability of each ant
            if s == (number_of_nodes+1) % last tabu, fill it with initial node
                tabu(:, end) = tabu(:, 1);
            else
                % determine the next node to be visited
                probability = zeros(1, number_of_nodes); % probability initialization, certain column in probability matrix = the probability in those city
                temp = colon(1, number_of_nodes); % matrix for unvisited node
                tabu1 = tabu(k, 1:end-1); % create clone of tabu
                tabu1(tabu1 == 0) = []; % remove visited node
                temp(tabu1) = []; % remove visited node
                for x = 1:size(temp, 2) % looping 1 until total nodes
                    probability(1, temp(x)) = tao(tabu(k, s-1), temp(x))^alpha * visibility(tabu(k, s-1), temp(x))^beta;
                end
                probability = probability./sum(probability);
                % choose the next node based on probability
                sorted = sort(probability);
                range = cumsum(sorted);
                randomResult = rand();
                founded1 = find(probability == sorted(max(find(range <= randomResult))+1));
                tabu(k, s) = founded1(floor((rand()*size(founded1, 2))+1)); % random when there're the same probabilities
            end
        end
    end
    
    % calculate the distance
    distance_of_ants_tour = zeros(number_of_ants, 1);
    for k=1:number_of_ants
        for x = 1:number_of_nodes
            distance_of_ants_tour(k, 1) = distance_of_ants_tour(k, 1) + round(distance(tabu(k, x), tabu(k, x+1)));
        end
    end
    
    % update best global solution (best of all iteration)
    founded = find(distance_of_ants_tour' == min(distance_of_ants_tour));
    if iteration == 1
        shortest_path = tabu(founded(1), :);
        shortest_distance_all_iteration(iteration) = min(distance_of_ants_tour);
    else
        if min(distance_of_ants_tour) < shortest_distance_all_iteration(iteration-1)
            shortest_path = tabu(founded(1), :);
            shortest_distance_all_iteration(iteration) = min(distance_of_ants_tour);
        else
            shortest_distance_all_iteration(iteration) = shortest_distance_all_iteration(iteration-1);
        end
    end
    
    % update best local solution (best of certain iteration)
    shortest_distance_each_iteration(iteration) = min(distance_of_ants_tour);
    
    % update delta tao of each iteration (global update)
    delta_tao = zeros(number_of_nodes, number_of_nodes);
    tao_addition_of_every_ant = Q ./ distance_of_ants_tour;
    for k = 1:number_of_ants
        for x = 1:(number_of_nodes)
            delta_tao(tabu(k, x), tabu(k, x+1)) = delta_tao(tabu(k, x), tabu(k, x+1)) + tao_addition_of_every_ant(k);
            delta_tao(tabu(k, x+1), tabu(k, x)) = delta_tao(tabu(k, x), tabu(k, x+1));
        end
    end

    % update tao
    % for k = 1:number_of_ants
    %     for x = 1:(number_of_nodes)
    %         tao(tabu(k, x), tabu(k, x+1)) = (tao(tabu(k, x), tabu(k, x+1)) .* rho) + delta_tao(tabu(k, x), tabu(k, x+1));
    %         tao(tabu(k, x+1), tabu(k, x)) = tao(tabu(k, x), tabu(k, x+1));
    %     end
    % end
    tao = (tao.*rho) + delta_tao;
    fprintf('iteration: %d/%d | shortest distance of all iteration: %d \n', iteration, MAX_ITERATION, shortest_distance_all_iteration(iteration));
    
    iteration = iteration + 1;
end

close
f = figure;
subplot(1, 2, 1);
plot([1:MAX_ITERATION]', shortest_distance_each_iteration(:, 1)');
title('Shortest Distance Each Iteration');
xlabel('Iterations');
ylabel('Shortest Distance');

subplot(1, 2, 2);
plot([1:MAX_ITERATION]', shortest_distance_all_iteration(:, 1)');
title('Shortest Distance All Iteration');
xlabel('Iterations');
ylabel('Shortest Distance');
drawnow
set(f, 'position', [200, 200, 700, 300]);
disp(shortest_path);
disp(shortest_distance_all_iteration(MAX_ITERATION));