function dYdt = SolveEquations(n,A,k,...
    nnodes,tau, GAMMA, max_growth,active_node_list,pause_relax)

dYdt=zeros(nnodes,1);

for i=1:nnodes
        
    if active_node_list(i) == 1
        
        if A(i)>0
            dIdt = n(i)*GAMMA(i,1);
            dRdt = n(i) * tau(i) * (1 - (n(i) / (k(i)))) * ((n(i) / A(i)) - 1);
            dYdt(i) = min([ (dRdt + dIdt) (n(i)*max_growth(i)) ]);
            
        elseif A(i) == 0
            dIdt = n(i)*GAMMA(i,1);
            dRdt = n(i) * tau(i) * (1 - (n(i) / (k(i)))) ;
            dYdt(i) = min([ (dRdt + dIdt) (n(i) * max_growth(i)) ]);
            
        elseif A(i) == -1
            dIdt = n(i) * GAMMA(i,1);
            dRdt = tau(i) * ( k(i) - n(i) ) ;
            dYdt(i) = dRdt + dIdt;
            
        else
        end
    end

    
% ------------ Note: Active / inactive nodes not used in v1.0.0 -----------
    if active_node_list(i)==0

        
        if pause_relax(i) == 1 % n --> k during dormancy
            if A(i) > 0
                dIdt = 0;
                dRdt = n(i) * tau(i) * (1 - (n(i) / (k(i)))) * ((n(i) / A(i)) -1);
                dYdt(i) = min([ (dRdt + dIdt) (n(i) * max_growth(i)) ]);
                
            elseif A(i) == 0
                dIdt = 0;
                dRdt = n(i) * tau(i) * (1 - (n(i) / (k(i)))) ;
                dYdt(i) = min([ (dRdt + dIdt) (n(i) * max_growth(i)) ]);
                
            elseif A(i) == -1
                dIdt = 0;
                dRdt = tau(i) * ( k(i) - n(i) ) ;
                dYdt(i) = dRdt + dIdt;
                
            elseif A(i) == -10 %(decay of Juveniles)
                
                dIdt = 0;
                dRdt = tau(i) * ( 0 - n(i) ) ;
                dYdt(i) = dRdt + dIdt;
            else
            end
            
        end
        
        if pause_relax(i) == 0  %n held constant during dormancy
            dYdt(i) = 0;
        end
        
    end
    
end

end



