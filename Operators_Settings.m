%% Operators configuration settings
% Pablo Caballero Garcés
% 30/03/15
function [ OpSettings ] = Operators_Settings(operators,s_o,belonging,NetSettings)

%% Number of operators
OpSettings.operators=operators; 

if size(s_o,2)==operators
    if size(belonging,2)==operators
        
        %% Network shares
        OpSettings.s_o=s_o; 
        
        %% Users per operator and ops_belong matrix
        ops_belongs=[];
        for i=1:operators
            N_k(i)=belonging(i);
            ops_belongs=[ops_belongs; i*ones(N_k(i),1)];
        end
        OpSettings.N_k=N_k; % vector of number of users per slice.
        OpSettings.ops_belongs=ops_belongs'; % user belonging vector, 
        % ops_belongs(i) = the idx of slice user i belongs to.
        
        %% Weights
        % According to equal weight sharing in SCPF,
        for i=1:NetSettings.users
            w_i(i)=s_o(ops_belongs(i))/N_k(ops_belongs(i));
        end
        OpSettings.w_i=w_i;
        %% Weights no sharing, deprecated.
        for i=1:NetSettings.users
            w_i_ss(i)=1/N_k(ops_belongs(i));
        end
        %OpSettings.w_i_ss=w_i_ss;        
        %% Checks
        if floor(100*sum(w_i)+0.000001)~=100
             error('Weight error')
        end
    else
        error('Belonging wrong input')
    end    
else
    error('Shares wrong input')
end

end
