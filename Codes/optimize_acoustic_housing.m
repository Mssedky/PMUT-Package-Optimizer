% Genetic algorithm optimization code for acoustic_cost.m 
% This code attempts to maximize the peak acoustic pressure of a PMUT and 
% logs the costs in a .txt file. You can validate and visualize the output
% designs in Solver_Acoustic.m and PMUT_Horn_GUI.m

% Code by Mostafa Sedky
% August 1, 2026


function S1 = optimize_acoustic_housing()  
    clc; clear all; close all;
    start_time = datetime;


    %% %% Process Optimization %% %%

    %% Definitions
    % Pi(Lambda)          anonymous function for evaluating fitness;
    %                          Pi.m should be in the working directory (in
    %                          this case we use acoustic_cost.m)
    % K,          1 x 1,  number of design strings to preserve and breed
    % TOL,        1 x 1,  cost function threshold to stop evolution
    % G,          1 x 1,  maximum number of generations
    % S,          1 x 1,  total number of design strings per generation
    % dv,         1 x 1,  number of design variables per string
    % PI,         G x S,  cost of sth design in the gth generation
    % Orig,       G x S,  indices of sorted strings before sorting
    %                     e.g. Orig(10, 1) = 34 means that the 1st ranked 
    %                          string in generation 10 was in position 34, 
    %                          visualize using familyTree.m
    % Lambda,     dv x S, array of most recent design strings
    % g,          1 x 1,  generation counter
    % PI_best,    1 x g,  minimum cost across strings and generations
    % PI_avg,     1 x g,  average cost across strings and generations
    % PI_par_avg, 1 x g,  average cost across strings and generations
    
    %% Givens
    K = 10;
    P = 10;
    TOL = 1e-7;
    G = 100;
    S = 30;
    dv = 11;
    
    horn_height_min = 1.0e-3;             % Height of acoustic housing (m)
    horn_height_max = 3.2e-3;             % Height of acoustic housing (m)
    r_top_min = 0.2e-3;                   % Top radius (previously 0.9e-3)
    r_top_max = 3.0e-3;                   % Top radius (previously 3.0e-3)
    p_min = 0.2;                       % Ellipsoid exponent p
    p_max = 10.0;                      % Ellipsoid exponent p
    q_min = 0.2;                       % Ellipsoid exponent q
    q_max = 10.0;                      % Ellipsoid exponent q
    r_min = 0.2;                       % Ellipsoid exponent r
    r_max = 10.0;                      % Ellipsoid exponent r
    A_theta_min = 0;                   % Azimuthal modulation amplitude
    A_theta_max = 0;                   % Azimuthal modulation amplitude
    m_theta_min = 0;                   % Azimuthal mode number
    m_theta_max = 0;                   % Azimuthal mode number
    A_x_min = 0;                       % X-direction modulation
    A_x_max = 0;                     % X-direction modulation
    A_y_min = 0;                       % Y-direction modulation
    A_y_max = 0;                     % Y-direction modulation
    A_z_min = 0;                       % Z-direction modulation
    A_z_max = 0;                     % Z-direction modulation
    m_z_min = 0;                       % Z-direction mode number
    m_z_max = 0;                       % Z-direction mode number
    
    % For construction of random genetic strings
    scale_factor = [horn_height_max - horn_height_min; ...
                r_top_max - r_top_min; ...
                p_max - p_min; ...
                q_max - q_min; ...
                r_max - r_min; ...
                A_theta_max - A_theta_min; ...
                m_theta_max - m_theta_min; ...
                A_x_max - A_x_min; ...
                A_y_max - A_y_min; ...
                A_z_max - A_z_min; ...
                m_z_max - m_z_min; ...
                ];


    offset = [horn_height_min; ...
                r_top_min; ...
                p_min; ...
                q_min; ...
                r_min; ...
                A_theta_min; ...
                m_theta_min; ...
                A_x_min; ...
                A_y_min; ...
                A_z_min; ...
                m_z_min; ...
                ];

    
    % Initialize
    PI = ones(G, S);
    Orig = ones(G, S);
    Lambda = rand(dv, S).*scale_factor + offset;

    % ===============================
    % Open log file for GA progress
    % ===============================
    log_filename = 'ga_convergence_log.txt';
    fid = fopen(log_filename, 'w');
    
    % Write header
    fprintf(fid, 'Generation\tMinCost\n');

    
    %% First generation
    g = 1;
    cost = ones(1,S);
    for i=1:S
        lambda_i = Lambda(:,i);
        % evaluate the fitness of each genetic string
        fprintf("Evaluating cost: gen=%d, string=%d\n\n", g, i);
        cost_i = acoustic_cost(lambda_i); 
        cost(i) = cost_i;
        fprintf("Parameters: L: %d, r_top: %d, p: %d, q: %d, r: %d, rest: %d, %d, %d, %d, %d, %d \n\n", lambda_i(1), lambda_i(2), lambda_i(3), lambda_i(4), lambda_i(5), lambda_i(6), lambda_i(7), lambda_i(8), lambda_i(9), lambda_i(10), lambda_i(11));
    end

    [new_cost, ind] = sort(cost); % order in terms of decreasing cost    
    PI(g, :) = new_cost;          % log the initial population costs
    Orig(g,:) = ind;              % log the indices before sorting
    Lambda = Lambda(:,ind);       % order in terms of decreasing cost
    
    % Store values for performance tracking
    PI_best = 1e10*ones(1,G);
    PI_avg = 1e10*ones(1,G); 
    
    top_performers = Lambda(:,1:4);
    top_costs = new_cost(1:4);
    
    % Update performance trackers
    PI_best(1) = min(new_cost);
    PI_avg(1) = mean(new_cost);
    MIN = min(new_cost);   
    
    % Log generation 1
    fprintf(fid, '%d\t%.8e\n', 1, PI_best(1));

    % Keep the sorted costs so we can reuse parent costs next generation
    prev_cost_sorted = new_cost;

    %% All later generations
    % (MIN > TOL) &&
    while  (g < G)
        g = g + 1;
        disp(strcat('Generation : ', num2str(g)))
        
        % Mating 
        parents = Lambda(:, 1:K);
        kids = zeros(dv, K);
        for p=1:2:K           % p = 1, 3, 5, 7,...      
            if mod(K, 2)
                disp('P is odd. Choose an even number of parents.')
                return
            end
            phi1 = rand(); 
            phi2 = rand();
            kids(:,p)   = phi1 * parents(:,p) + (1 - phi1) * parents(:,p+1);
            kids(:,p+1) = phi2 * parents(:,p) + (1 - phi2) * parents(:,p+1);
        end
        new_strings = rand(dv, S-2*K).*scale_factor + offset;
        
        % Update Lambda
        Lambda = [parents, kids, new_strings]; % concatenate horizontally
        
        % Evaluate fitness of new population        
        cost = ones(1,S);
        cost(1:K) = prev_cost_sorted(1:K);
        % compute only for kids and new random strings
        for i=K+1:S
            lambda_i = Lambda(:,i);
            fprintf("Evaluating cost: gen=%d, string=%d\n\n", g, i);
            cost_i = acoustic_cost(lambda_i); 
            cost(i) = cost_i;
            fprintf("Parameters: L: %d, r_top: %d, p: %d, q: %d, r: %d, rest: %d, %d, %d, %d, %d, %d \n\n", lambda_i(1), lambda_i(2), lambda_i(3), lambda_i(4), lambda_i(5), lambda_i(6), lambda_i(7), lambda_i(8), lambda_i(9), lambda_i(10), lambda_i(11));
        end
        
        [new_cost, ind] = sort(cost);     
        PI(g, :) = new_cost;        
        Orig(g,:) = ind; 
        Lambda = Lambda(:,ind);
        % disp(Lambda(1,1:6))
    
        % Update performance trackers
        PI_best(g) = min(new_cost);
        PI_avg(g) = mean(new_cost);
        
        if min(new_cost) < MIN
            MIN = min(new_cost);
        end

         % Save sorted costs for next generation reuse
        prev_cost_sorted = new_cost;

        top_costs = new_cost(1:4);
        top_performers = Lambda(:,1:4);

        % Log to file
        fprintf(fid, '%d\t%.8e\n', g, PI_best(g));
        fprintf("Generation %d best cost = %.6e\n \n\n", g, PI_best(g));
    end

    % Close log file
    fclose(fid);
    
    end_time = datetime;
    sim_time = between(start_time, end_time);
    disp(sim_time)

    %% %% DELIVERABLES %% %%
  
    % Plot the cost function over time, demonstrate convergence
  
    figure(1)
    semilogy(2:g, PI_best(2:g), 'LineWidth', 2);
    % hold on
    % semilogy(2:g, PI_avg(2:g), 'LineWidth', 2);
    xlabel('Generations', 'Interpreter', 'latex', 'FontSize', 20);
    ylabel('Cost', 'Interpreter', 'latex', 'FontSize', 20);
    title('pMUT Acoustic Housing: Convergence of Cost Function', 'Interpreter', 'latex', 'FontSize', 20);
    % legend('Best', 'Overall Mean', 'Parent Mean', 'Interpreter', 'latex', 'FontSize', 15, 'location', 'east');
    
    % Report 4 best-performing designs
    S1 = top_performers(:,1);
    disp('Top Performing Strings')
    disp('S1')
    disp(top_performers(1,1))
    disp(top_performers(2,1))
    disp(top_performers(3,1))
    disp(top_performers(4,1))
    disp(top_performers(5,1))
    disp(top_performers(6,1))
    disp(top_performers(7,1))
    disp(top_performers(8,1))
    disp(top_performers(9,1))
    disp(top_performers(10,1))
    disp(top_performers(11,1))
    disp('S2')
    disp(top_performers(1,2))
    disp(top_performers(2,2))
    disp(top_performers(3,2))
    disp(top_performers(4,2))
    disp(top_performers(5,2))
    disp(top_performers(6,2))
    disp(top_performers(7,2))
    disp(top_performers(8,2))
    disp(top_performers(9,2))
    disp(top_performers(10,2))
    disp(top_performers(11,2))
    disp('S3')
    disp(top_performers(1,3))
    disp(top_performers(2,3))
    disp(top_performers(3,3))
    disp(top_performers(4,3))
    disp(top_performers(5,3))
    disp(top_performers(6,3))
    disp(top_performers(7,3))
    disp(top_performers(8,3))
    disp(top_performers(9,3))
    disp(top_performers(10,3))
    disp(top_performers(11,3))

    disp('S4')
    disp(top_performers(1,4))
    disp(top_performers(2,4))
    disp(top_performers(3,4))
    disp(top_performers(4,4))
    disp(top_performers(5,4))
    disp(top_performers(6,4))
    disp(top_performers(7,4))
    disp(top_performers(8,4))
    disp(top_performers(9,4))
    disp(top_performers(10,4))
    disp(top_performers(11,4))
    
    disp('Costs of Top Performing Strings')
    disp('Pi1 | Pi2 | Pi3 | Pi4')
    disp(top_costs)
  
end