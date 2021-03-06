%% Fair belonging
% Pablo Caballero Garc??s
% 30/03/15
function [ belonging ] = get_belonging(s_o,users,operators)
% Get the number of users belonging to each slice.
% Parameters:
% s_o: vector of share per slice
% users: total number of users
% operators: number of operators
% Returns:
% belonging: vector, belonging(i) is the user number under slice i.
for op=1:operators-1
    belonging(op)=round(s_o(op)*users);
end
belonging(operators)=users-sum(belonging);
end