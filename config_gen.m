function savefile = config_gen(prefix,n_LHS, config_dir, plot_network, plot_timeseries)

fprintf('----------------  O S I R I S  v1.0.0  ----------------');

fprintf('\n\n\nCONFIGURATION FILE PRODUCTION:');
fprintf('\n\tReading input files...');

cd(config_dir);
%Run conditions
basic = csvread('basic.csv');

%Node properties
FileName = 'node_names.csv';
fileID = fopen(FileName,'r');
formatSpec = '%s';
N_names = textscan(fileID, formatSpec, 'ReturnOnError', false); fclose(fileID);
tau = csvread('tau.csv');
A = csvread('allee_A.csv');
initial_n = csvread('initial_n.csv');
max_growth = csvread('max_growth.csv');
max_growth_errors = csvread('max_growth_errors.csv');
prey_switch = csvread('prey_switch.csv');

%Network properties
adjacency_matrix = csvread('adjacency_matrix.csv');

%Interactions: n-->n source, target and interaction type
y_int_source_target = csvread('y_int_source_target.csv');
y_int = csvread('y_int.csv');%4 parameters for function n->n
y_int_error = csvread('y_int_error.csv');%

%External forcing
FileName = 'forcing_names.csv'; 
fileID = fopen(FileName,'r');
formatSpec = '%s';
F_names = textscan(fileID, formatSpec, 'ReturnOnError', false); fclose(fileID);
forcing_names = vertcat(F_names{:});
forcing = csvread('forcing.csv');
forcing_noise_magnitude = csvread('forcing_noise_magnitude.csv'); %stdev of noise for each forcing
forcing_stdev = csvread('forcing_stdev.csv'); %parameters for linear trend in stdev of each forcing

%Interactions: F-->n source, target and interaction type
y_ext_source_target = csvread('y_ext_source_target.csv');
y_ext = csvread('y_ext.csv');%4 parameters for function F-->n
y_ext_error = csvread('y_ext_error.csv');%errors on each parameter

%Interactions: F-->k source, target and interaction type
phi_source_target = csvread('phi_source_target.csv');
phi_Fm = csvread('phi_Fm.csv');%4 parameters for function F-->k
phi_Fm_error = csvread('phi_Fm_error.csv');%errors on each parameter
chi_Fk = csvread('chi_Fk.csv'); %additive and multiplicative coefficients
phi_Fm_ref = csvread('phi_Fm_ref.csv');%redundant: use source/target instead

%Run conditions
nthdatapoint = basic(9);
duration = basic(5);
output_resolution = basic(6);
n_iterations = duration / output_resolution; %number of reporting points

%Perturbation ('crisis') events
include_crisis = basic(10);
nodes_affected = csvread('nodes_affected.csv');
crisis_times = csvread('crisis_times.csv');
crisis_values = csvread('crisis_values.csv');
crisis_value_definition = basic(11);

%Special case of Oxygen node (not specified/used in present version)
oxygen_nodes = 0; %csvread('oxygen_nodes.csv');
oxygen_thresholds = 0; %csvread('oxygen_thresholds.csv');
include_oxygen_thresholds = 0; %basic(12);
oxygen_state_node = 0; %basic(13);

cd ..; %back to base directory after reading config files

fprintf(' done.\n');

%---------------- Define parameters ---------------------------------------
%Network
nnodes = basic(1);
node_names = vertcat(N_names{:});
n_F = basic(7);
n_interactionsFk = basic(4);
n_interactionsFn = basic(3);
n_interactionsN = basic(2);
n_interactions = sum(sum(adjacency_matrix));

%Interactions: n-->n source, target and interaction type
y_int_source = y_int_source_target(:,1); %source_nodes
y_int_target = y_int_source_target(:,2); %N_target_nodes
y_int_type   = y_int_source_target(:,3); %N_interact_type


%EXTERNAL FORCING
t_F_import = 0:length(forcing); %real time yrs
dt = duration / n_iterations;
t_F = linspace(0, duration, n_iterations);
F_t = zeros(numel(t_F),n_F);

for i=1:n_F
    F_t(:,i) = interp1(t_F_import(1:duration+1), forcing(1:duration+1,i), t_F);
end


%Interactions: F-->n source, target and interaction type
y_ext_source = y_ext_source_target(:,1);
y_ext_target = y_ext_source_target(:,2);
y_ext_type   = y_ext_source_target(:,3);

%Interactions: F-->k source, target and interaction type
phi_source = phi_source_target(:,1);
phi_target = phi_source_target(:,2);
phi_type   = phi_source_target(:,3);

PHI=zeros(nnodes,1);

%Calculate size of predation groups, for prey-switching
[ng, prey_group_members] = prey_groups(y_int_source_target, y_int, oxygen_state_node, nnodes);


%---------------- LHS process ---------------------------------------------
%create pristine copies of values
copy_phi_Fm = phi_Fm;
copy_y_ext = y_ext;
copy_y_int = y_int;

%number and position of variables to sample
phi_Fm_pos = find(phi_Fm_error > 0);
phi_FM_n = numel(phi_Fm_pos); %phi_Fm

y_ext_pos = find(y_ext_error > 0);
y_ext_n = numel(y_ext_pos); %'y_ext(x)'

y_int_pos = find(y_int_error > 0);
y_int_n = numel(y_int_pos); %'y_int(x)'
total = phi_FM_n + y_ext_n + y_int_n; %total number of LHS dimensions needed

% LHS
L = lhsdesign(total, n_LHS);

%do selective modification and file creation
for j=1:n_LHS
    
    output_dir = 'Output_files';
    run_id=j;
    
    %set pristine copies
    phi_Fm=copy_phi_Fm;
    y_ext=copy_y_ext;
    y_int=copy_y_int;
     
    %add noise to forcings
    if j > 1
        for i = 1:n_F
            %take annual values
            annual_values = forcing(1:duration + 1, i);
            %create noise for each value
            rnd_z = normrnd(0, forcing_noise_magnitude(i), [length(annual_values), 1]);
            %add noise to original values
            noisy_annual_values = annual_values + rnd_z;
            %create new interpolated version of the randomized annual values
            original_time = t_F_import(1:duration+1)';
            F_t(:,i) = interp1(original_time, noisy_annual_values, t_F');
            
        end
    end
    
    %Phi_Fm - calculate variants
    for i = 1:phi_FM_n
        upper = copy_phi_Fm(phi_Fm_pos(i)) + phi_Fm_error(phi_Fm_pos(i));
        lower = copy_phi_Fm(phi_Fm_pos(i)) - phi_Fm_error(phi_Fm_pos(i));
        a = lower + (upper-lower) * L(i,j);
        phi_Fm(phi_Fm_pos(i)) = a; %change relevant element in copy
        
        if j == 1
            phi_Fm(phi_Fm_pos(i)) = copy_phi_Fm(phi_Fm_pos(i));
        end
        
    end
    
    %y_ext - calculate variants
    for i = 1:y_ext_n
        upper = copy_y_ext(y_ext_pos(i)) + y_ext_error(y_ext_pos(i));
        lower = copy_y_ext(y_ext_pos(i)) - y_ext_error(y_ext_pos(i));
        a = lower + (upper-lower) * L(i + phi_FM_n,j);
        y_ext(y_ext_pos(i)) = a; %change relevant element in copy
        if j == 1
            y_ext(y_ext_pos(i)) = copy_y_ext(y_ext_pos(i));
        end
    end
    
    %y_int - calculate variants
    for i = 1:y_int_n
        upper = copy_y_int(y_int_pos(i)) + y_int_error(y_int_pos(i));
        lower = copy_y_int(y_int_pos(i)) - y_int_error(y_int_pos(i));
        a = lower + (upper-lower) * L(i + phi_FM_n + y_ext_n, j);
        y_int(y_int_pos(i)) = a; %change relevant element in copy
        if j == 1
            y_int(y_int_pos(i)) = copy_y_int(y_int_pos(i));
        end
    end
    
    %save workspace in file
    savefile = [prefix '_' int2str(j) '.mat'];
    fprintf('\tConfiguration file produced: %s\n', savefile);
    save(savefile);

    
    %for batch runs:
    %    savefile = [d_name '/config' prefix '_' int2str(j) '.mat'];
    %    save(savefile, '-regexp', '^(?!(L)$).');
    
end

end
