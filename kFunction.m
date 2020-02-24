function[y] = kFunction(tol_func_type,a,b,c,d,F_mean,F_mean_ref,F_sig,F_sig_ref)

%Forcing
forcing_mean = [F_mean_ref, F_mean]; %first value is reference point
forcing_sigma = [F_sig_ref, F_sig];

%Tolerance (remains fixed over time - potential to change this in future)
T_a = [a, a]; %(in normalized units)
T_b = [b, b];
T_c = [c, c];
T_d = [d, d];
G_T = 0;

sum_tol = zeros(2,1);

for i=1:2
    
    %create forcing distribution
    lower_F = forcing_mean(i) - (4 * forcing_sigma(i));
    upper_F = forcing_mean(i) + (4 * forcing_sigma(i));
    x = linspace(lower_F, upper_F,100); %range of relative forcing values (for T and F)
    G_F = ( exp(-((x - forcing_mean(i)) ./ forcing_sigma(i)).^2));
    
    %create tolerance response
    if tol_func_type == 1 %Linear
        G_T = T_a(i) .* x + T_b(i);
        
    elseif tol_func_type==2 %Logistic
        G_T = T_a(i) ./ (1 + exp(-T_b(i) .* (T_c(i) - x))) + T_d(i);
        
    elseif tol_func_type == 3 %Polynomial
        G_T = T_a(i) + T_b(i).* x + (T_c(i) .* x.^2) + (T_d(i) .* x.^3);
        
    elseif tol_func_type == 4 %Gaussian 
        G_T = ( T_a(i) .* exp(-(((x - T_b(i)).^2) ./ (2 * T_c(i).^2)))) + T_d(i);    
        
    elseif tol_func_type == 5 %Slope-only, with y=1 when x=1
        G_T = T_a(i) .* x + (1 - T_a(i));
        
    elseif tol_func_type == 6 %Slope-only, with y=0 when x=1
        G_T = T_a(i) .* x + (-T_a(i));
        
    elseif tol_func_type == 7 %Slope-only, with y=0 when x=0
        G_T = T_a(i) .* x;
        
    elseif tol_func_type == 8 %Slope-only, with y=1 when x=1
        G_T = T_a(i) .* x + (1 - T_a(i));         
        
    end
    
    %multiply & integrate
    G_FG_T = G_F .* G_T;
    sum_tol(i) = sum(G_FG_T);
     
end

y = sum_tol(2) / sum_tol(1); %normalize to standard condition

end

