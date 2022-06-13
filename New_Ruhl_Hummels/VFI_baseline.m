function [V] = VFI_baseline(V_NX, V_EX, V_CX, rho_EE, rho_CC, trans_joint, R)
% Perform value function iteration for baseline sunk cost model
iterate = 0;
error = 1;
while error > 1e-5 && iterate < 2000
    V_NX(:, 8) = V_NX(:, 5) + R * (trans_joint * V_NX(:, 7));
    V_NX(:, 9) = V_NX(:, 6) + R * (trans_joint * V_EX(:, 7));
    V_EX(:, 8) = V_EX(:, 5) + R * (trans_joint * V_NX(:, 7));
    V_EX(:, 9) = V_EX(:, 6) + R * (rho_EE * (trans_joint * V_EX(:, 7)) + (1 - rho_EE) * (trans_joint * V_CX(:, 7)));
    V_CX(:, 8) = V_CX(:, 5) + R * (trans_joint * V_NX(:, 7));
    V_CX(:, 9) = V_CX(:, 6) + R * (rho_CC * (trans_joint * V_CX(:, 7)) + (1 - rho_CC) * (trans_joint * V_EX(:, 7)));
    % Maximize value function
    V_NX(:,11) = max(V_NX(:, 8), V_NX(:, 9));
    V_EX(:,11) = max(V_EX(:, 8), V_EX(:, 9));
    V_CX(:,11) = max(V_CX(:, 8), V_EX(:, 9));
    % Get error
    a = max(abs(V_NX(:, 7) - V_NX(:, 11)));
    b = max(abs(V_EX(:, 7) - V_EX(:, 11)));
    c = max(abs(V_CX(:, 7) - V_CX(:, 11)));
    error = max([a b c]);
    % update value function
    V_NX(:, 7) = V_NX(:, 11); 
    V_EX(:, 7) = V_EX(:, 11); 
    V_CX(:, 7) = V_CX(:, 11); 
    iterate = iterate + 1;
end

% Combine the three matrix
V = [V_NX;V_EX;V_CX];
loops = length(V);
% Policy function
for j = 1:loops
    if V(j, 8) > V(j, 9)
        V(j, 10) = 0;
    else
        V(j, 10) = 1;
    end
end
end

