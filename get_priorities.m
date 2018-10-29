%% Priorities assignment
% Pablo Caballero Garces
% 08/06/16

function [ phi ] = get_priorities(OpSettings,levels,users)
phi=[];
if levels==1
    for op=1:OpSettings.operators
        add=ones(1,OpSettings.N_k(op));
        %add=OpSettings.s_o(op).*add./sum(add);
        add=add./sum(add);
        phi=[phi add];
    end
else
    if levels==2
        for op=1:OpSettings.operators
            add=[ones(1,floor(OpSettings.N_k(op)/2))...
                2*ones(1,-floor(OpSettings.N_k(op)/2)+OpSettings.N_k(op))];
            %add=OpSettings.s_o(op).*add./sum(add);
            add=add./sum(add);
            phi=[phi add];
        end
    elseif levels==4
        for op=1:OpSettings.operators
            steps=[floor(OpSettings.N_k(op)/4) floor(OpSettings.N_k(op)/4) ...
                floor(OpSettings.N_k(op)/4) OpSettings.N_k(op)-3*floor(OpSettings.N_k(op)/4)];
            add=[ones(1,steps(1)) 2*ones(1,steps(2)) 3*ones(1,steps(3))...
                4*ones(1,steps(4))];
            %add=OpSettings.s_o(op).*add./sum(add);
            add=add./sum(add);
            phi=[phi add];
        end
    else
        error('Wrong level')
    end
end
end