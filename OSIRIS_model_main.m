function [ndata_output] = OSIRIS_model_main(infile, LHS_i, n_LHS)

fprintf('\n\nMODEL RUN:');

if n_LHS > 1
    fprintf(' %d/%d', LHS_i, n_LHS)
end

duration = 0;

input_file = fullfile(pwd, infile);

fprintf('\n\tReading configuration file... ');
load(input_file);
fprintf('done.');

%Initialize 
k0 = ones(nnodes, 1); 
k0_phi_add = ones(nnodes, n_F); 
k1_phi_add = ones(nnodes, n_F);

gamma0Fn = zeros(nnodes, n_F);
gamma1Fn = zeros(nnodes, n_F);
gamma0nn = zeros(nnodes, nnodes);
gamma1nn = zeros(nnodes, nnodes);
GAMMA = zeros(nnodes, 1);

PHI = ones(nnodes, 1);
data_point_counter = 0;
F_t_int = zeros(1, n_F);

%Integration options
options = odeset('AbsTol', 1e-10);
integration_interval = [0, output_resolution];

%Data storage
n_points = int32(n_iterations / nthdatapoint) + 1;
n_data = zeros(n_points, nnodes);
t_data = zeros(n_points, 1);

%store initial conditions at t=0
n_data(1, 1:nnodes) = initial_n'; 
t_data(1, 1) = 0;


fprintf('\n\tRunning model...');
f = waitbar(0,'1','Name','Model running...',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
cancelled = false();

% main calculation loop
tic();

for iter = 1:n_iterations
    
    if getappdata(f, 'canceling')
        cancelled = true();
        fprintf('\n\n>>> CANCELLED <<<\n\n');
        break
        
    end
    
    data_point_counter = data_point_counter + 1;
    
    time_interval = [((iter - 1) * output_resolution) (iter * output_resolution)];
    tx = iter * output_resolution;
    
    for i = 1:n_F
        if iter <= 1 || iter >= (n_iterations - 1)
            F_t_int(1, i) = F_t(iter, i);
        else
            %locally interpolate forcing values
            F_t_int(1, i) = interp1 ( t_F(iter-1:iter + 1) , F_t(iter-1:iter + 1, i) , tx); 
        end
    end
   
    %calculate forcing stdev (inputs are relative stdev values)
    for i = 1:n_F
        F_stdev_int(1, i) = forcing_stdev(1, i) + forcing_stdev(2, i) * tx;
    end
    
    forcing_stdev_ref = forcing_stdev(1,:);
    
    additive_k = 0; 
    multiplicative_k = 0; 
    additive_g_nn = 0; 
    additive_g_Fn = 0;
    
    for j=1:nnodes
        
        %EXTERNAL EFFECT :  F -> k
        x = find(phi_target == j);
        nx = length(x);
        
        additive_k = 0; 
        multiplicative_k = 0; 
        additive_g_nn = 0; 
        additive_g_Fn = 0;
        
        k1 = zeros(n_F, 1);
        
        if nx >= 1
            %additive part
            for p = 1:nx
                %Forcing parameters
                F = F_t_int(1, phi_source(x(p))); %forcing value
                F_ref = 1;
                Fstdev = F_stdev_int(1, phi_source(x(p))); %forcing stdev value
                Fstdev_ref = forcing_stdev_ref(phi_source(x(p))); %forcing stdev value
                %note: j==phi_target(x(p))
                
                %Tolerance parameters
                row = (j - 1) * n_F + phi_source(x(p));
                a = phi_Fm(row, 1);
                b = phi_Fm(row, 2);
                c = phi_Fm(row, 3);
                d = phi_Fm(row, 4);
                
                k1_phi_add(j, phi_source(x(p))) = kFunction(phi_type((x(p))),a,b,c,d, F,F_ref,Fstdev, Fstdev_ref);
                
                additive_k = additive_k + ((chi_Fk(row, phi_source(x(p)))) *...
                    (k1_phi_add(j, phi_source(x(p))) - k0_phi_add(j, phi_source(x(p)))));
            end
            
            %multiplicative part
            if nx > 1 %needs to be >1 factors in order to multiply their effects

                permutations = nchoosek(unique(phi_source(x)), 2);
                [n_permutations, ~] = size(permutations);
                
                for q = 1:n_permutations
                    chi_temp = chi_Fk((j - 1) * n_F + permutations(q, 1), permutations(q, 2));
                    k1_phi_mult_F1 = k1_phi_add(j, permutations(q, 1));
                    k1_phi_mult_F2 = k1_phi_add(j, permutations(q, 2));
                    multiplicative_k = multiplicative_k + (chi_temp*...
                        (((k1_phi_mult_F1 / k0_phi_add(j, permutations(q, 1)))...
                        * (k1_phi_mult_F2 / k0_phi_add(j, permutations(q, 2)))) - 1));
                end
            end
        end
        
        x_k = k0(j) + additive_k + (multiplicative_k * k0(j));
        PHI(j) = max(0, x_k);
                
        %EXTERNAL EFFECT :  F --> n
        x = find(y_ext_target == j);
        nx = length(x);
        
        %additive part (no multiplicative part for F-->n)
        if nx >= 1
            for p = 1:nx    %for each relevant forcing
                F = F_t_int(1, y_ext_source(x(p))); %forcing value
                row = (j - 1) * n_F + y_ext_source(x(p)); %empirical parameters
                a = y_ext(row, 1);
                b = y_ext(row, 2);
                c = y_ext(row, 3);
                d = y_ext(row, 4);

                y_x = InteractionFunction(y_ext_type(x(p)),a,b,c,d,F);

                if A(j) > 0 %note: y_ext_target(x(p))==j
                    gamma1Fn(j, y_ext_source(x(p))) = y_x;
                else
                    gamma1Fn(j, y_ext_source(x(p))) = y_x;
                end
                additive_g_Fn = additive_g_Fn + (1 * (gamma1Fn(j, y_ext_source(x(p)))...
                    - gamma0Fn(j, y_ext_source(x(p)))));
            end
        end
        
        
        %INTERNAL EFFECT : n --> n
        
        intercept_correction = intercept_adjustment(nnodes, ng,prey_switch, initial_n, prey_group_members, y_int);
        
        %find nodes targeting j
        x = find(y_int_target == j); %note: y_int_target(x(p))==j
        nx = length(x); 
        
        %additive part (no multiplicative part for n->n)
        if nx >= 1
            for p = 1:nx  %for each relevant forcing
                n_i = initial_n(int32((y_int_source(x(p))))); %value of source node
                row = (j - 1) * 40 + y_int_source(x(p));
                a = y_int(row,1);
                b = y_int(row,2);
                c = y_int(row,3);
                d = y_int(row,4);
                
                %caluculate interaction
                correction_value = intercept_correction((j - 1) * nnodes + y_int_source(x(p)), 3) * -a;
                y_x = InteractionFunction(y_int_type(x(p)),a,b,c,d,n_i) + correction_value; 
                
                gamma1nn(j, y_int_source(x(p))) = y_x; 
                
                additive_g_nn=additive_g_nn + (1 * (gamma1nn(j, y_int_source(x(p)))...
                    - gamma0nn(j, y_int_source(x(p))) )); 
            end
        end
        
        %combine F-->n and n-->n effects
        GAMMA(j) = GAMMA(j) + additive_g_nn + additive_g_Fn;
        
    end %j
    
    
    %Perturbation ('crisis') events
    if include_crisis == 1
        if sum(ismember(crisis_times, iter * output_resolution)) > 0
            if crisis_value_definition == 1
                for cr = 1:length(nodes_affected)
                    if iter * output_resolution == crisis_times(cr)
                        initial_n(nodes_affected) = crisis_values(cr);
                    end
                end
            end
            if crisis_value_definition == 2
                
                for cr = 1:length(nodes_affected)
                    if iter * output_resolution == crisis_times(cr)
                        initial_n(nodes_affected(cr)) = initial_n(nodes_affected(cr)) * crisis_values(cr);
                    end
                end
            end
        end
    end

    %un-used in this version, included to allow compatability with newer version of 'SolveEquations.m'
    active_node_list = ones(nnodes,1);
    pause_relax = ones(nnodes,1); 
    
    
    % Numerical integration
    [~, n] = ode45(@(t,n) SolveEquations(n,A,PHI,nnodes,tau, GAMMA, max_growth,...
        active_node_list, pause_relax), integration_interval,initial_n,options);
    
  
    if data_point_counter == nthdatapoint
        n_data(int32(iter / nthdatapoint) + 1, 1:nnodes) = n(end, :);
        t_data(int32(iter / nthdatapoint) + 1, 1) = time_interval(2);
        data_point_counter = 0;
        
        % Update waitbar
        so_far = iter / n_iterations;
        waitbar(so_far,f,sprintf('Progress:')) 
        
    end
    
    %arbitrary test for unstable solutions
    find_unstable = n(end,:) > 20; %obviously extreme population
    if sum(find_unstable) > 0
        n(end,:) = 1; %avoid explosion crash
        disp('>>>> Unstable solution found, check results <<<<'); %Alert
    end
       
    % reset as necessary
    initial_n = n(end,:);
    gamma0Fn = gamma1Fn;
    gamma0nn = gamma1nn;
    k0_phi_add = k1_phi_add;
    k0 = PHI;
 
end %iter

delete(f);

if ~cancelled
    fprintf('done. (model run time = %.1fs)\n', toc());
    
    %save data to file
    fprintf('\tSaving output files...');
    out_filename = strcat(output_dir, prefix, '_', int2str(run_id), '.csv');
    out_file = fullfile(pwd, output_dir, char(out_filename));
    dlmwrite(out_file,n_data,'precision','%.6f');
    
    fprintf('done.\n')
    ndata_output = n_data;
end


%% Plot timeseries

if plot_timeseries && ~cancelled
    
    if LHS_i == n_LHS
        fprintf('\n\nFIGURES:');
        fprintf('\n\tProducing timeseries figures...');
        
        for figs=1:2
            
            figure (figs);
            
            for k=1:4
                nnlist = [ 4 4 4 5  ]; % Note: this would need adjusting for new models
                nlist=   [  1  2  3  4  0   0  0  ;...
                    5  6  7  8  0   0  0  ;...
                    9  10 11 12 0   0  0  ;...
                    13 14 15 16 17  0  0  ];
                
                subplot(2,2,k);
                plot(t_data(:,1),n_data(:,nlist(k,1:nnlist(k))),'-');
                
                legend(node_names(nlist(k,1:nnlist(k))),'location','Northeast','FontSize',12);
                legend BOXOFF;
                set(gcf, 'units','normalized','outerposition',[0.0 0.0 1.2 1.2]);
                set(gca,'fontsize',16);
                
                if figs == 2
                    axis([0 duration 0 2]);
                end
                
            end
        end
        
        fprintf('done.')
    end
    
end


%% Plot networks

if plot_network && ~cancelled
    
    if LHS_i == n_LHS
        
        fprintf('\n\tProducing networks figures...');
        
        s = [y_int_source(1:n_interactionsN)' y_ext_source(1:n_interactionsFn)' + nnodes...
            phi_source(1:n_interactionsFk)' + nnodes];
        t = [y_int_target(1:n_interactionsN)' y_ext_target(1:n_interactionsFn)'...
            phi_target(1:n_interactionsFk)' ];
        
        s1 = y_int_source(1:n_interactionsN)';
        t1 = y_int_target(1:n_interactionsN)';
        
        s2 = y_ext_source(1:n_interactionsFn)' + nnodes;
        t2 = y_ext_target(1:n_interactionsFn)';
        
        s3 = phi_source(1:n_interactionsFk)' + nnodes;
        t3 = phi_target(1:n_interactionsFk)' ;
        
        names = [node_names;forcing_names];
        p = nnodes+n_F;
        
        figure(3);
        weights = ones(1, numel(t1));
        G = digraph(s1, t1,weights, p);
        h = plot(G,'NodeLabel',names,'Layout','circle','EdgeColor','k');
        h.NodeFontSize = 16;
        G.Nodes.Size = indegree(G);
        axis off; 
        axis tight; 
        pbaspect([1 1 1]);
        highlight(h,1:nnodes); 
        highlight(h,1:nnodes,'NodeColor','g');
        highlight(h,nnodes+1:nnodes+n_F);
        highlight(h,nnodes+1:nnodes+n_F,'NodeColor','b');
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        
        figure(4)
        weights = ones(numel(t2) + numel(t3),1);
        G = digraph([s2,s3],[t2,t3],1,p);
        h = plot(G,'NodeLabel',names,'Layout','circle');
        h.NodeFontSize = 16;
        
        highlight(h, 1:nnodes); 
        highlight(h, 1:nnodes,'NodeColor','g');
        highlight(h, nnodes + 1:nnodes + n_F);
        highlight(h, nnodes + 1:nnodes + n_F,'NodeColor','b');
        highlight(h,[s2],[t2],'EdgeColor','k');
        axis off; 
        axis tight; 
        pbaspect([1 1 1]);
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        
        fprintf('done.\n')
    end
    
end

if LHS_i == n_LHS
    fprintf('\n\nFINISHED.\n\n');
    fprintf('-------------------------------------------------------\n\n\n');
end

end
