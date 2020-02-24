function [ng, prey_group_members] = prey_groups(y_int_source_target, y_int, oxygen_state_node, nnodes)

%%% Clip off zeros
del = y_int_source_target == 0;
y_int_source_target(del) = [];
del = int32(length(y_int_source_target) / 3);
y_int_source_target = reshape(y_int_source_target,del, 3);

[np, ~] = size(y_int_source_target);

% Produce prey table
negs = (find(y_int(:,1) > 0));
n1 = length(negs);
p1_table = zeros(n1, 3);

a_source = 0;
a_targ = 0;
a_val = 0;
counter = 0;

prey_group_members = zeros(nnodes, nnodes);

for i = 1:np
    a_source = y_int_source_target(i, 1); 
    a_targ = y_int_source_target(i, 2);
    a_val = y_int(((a_targ) - 1) * 40 + a_source,1);
    
    if(a_val > 0)
        if a_source ~= oxygen_state_node
            counter = counter + 1;
            p1_table(counter, 1) = a_source;
            p1_table(counter, 2) = a_targ;
            p1_table(counter, 3) = a_val;
        end
    end
end

if n1 > counter
    del = p1_table == 0;
    p1_table(del) = [];
    p1_table = reshape(p1_table, (length(p1_table) / 3), 3);
end

p2_table = sortrows(p1_table, 2);

%%% Find group sizes
ng = zeros(nnodes, 1);

for i = 1:nnodes
    pos = find(p2_table(:, 2) == i);
    
    if isempty(pos)
        ng(i) = 0;
    else
        ng(i) = numel(pos); %fill in count array & prey_group_members matrix
        members = p2_table(pos, 1);
        prey_group_members(i, 1:ng(i)) = members';
    end    
end
end

