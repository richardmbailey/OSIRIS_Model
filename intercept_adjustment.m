function [intercept_correction] = intercept_adjustment(nnodes, ng,prey_switch, initial_n, prey_group_members, y_int)

intercept_correction=zeros(nnodes^2,3); %default value zero if no prey switching

for i = 1:nnodes
    for j = 1:nnodes
        intercept_correction((i - 1) * nnodes + j, 1) = i;
        intercept_correction((i - 1) * nnodes + j, 2) = j;
    end
end

for j = 1:nnodes
    
    Y0 = zeros(ng(j), 1); %intercepts under baseline equilibrium case (all n=1)
    Ystar = zeros(ng(j), 1);
    
    if prey_switch(j) > 0
        total_prey = 0;
        sum_Ystar = 0;
        sum_Y0 = 0;
        Yi = zeros(ng(j), 1);
        dYi = zeros(ng(j), 1);
        
        for tp = 1:ng(j) %total prey
            total_prey = total_prey + initial_n(prey_group_members(j, tp));
        end
        
        if total_prey > 0
            
            for tp = 1:ng(j) %sums for calculation of intercept correction
                current_prey_node = prey_group_members(j, tp);
                current_prey_state = initial_n(current_prey_node);
                
                Y0(tp,1) = -y_int(((j - 1) * 40) + prey_group_members(j, tp), 1); %(negative of slope)
                
                sum_Y0 = sum_Y0 + Y0(tp, 1);
                
                Ystar(tp,1) = Y0(tp, 1) - Y0(tp, 1) * ((current_prey_state/total_prey) - 1);
                sum_Ystar = sum_Ystar + Ystar(tp, 1);
                
            end
            
            
            for tp = 1:ng(j) %calculate intercept corrections
                
                Yi(tp, 1) = ((Ystar(tp, 1) / sum_Ystar)) * sum_Y0;
                dYi(tp,1) = (Yi(tp, 1) / Y0(tp, 1)) - 1;
                
            end
            
            for tp = 1:ng(j) %calculate intercept correction
                current_prey_node = prey_group_members(j, tp);
                intercept_correction((j - 1) * nnodes + current_prey_node, 3) = dYi(tp, 1) * prey_switch(j); %predator benefits
                intercept_correction((current_prey_node - 1) * nnodes + j, 3) = dYi(tp, 1) * prey_switch(j); %prey costs
                
            end
            
        end
    end
    
end

end