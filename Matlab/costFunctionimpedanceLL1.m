% Equivalent Circuit Model Function
function cost = costFunctionimpedanceLL1(x, C, w, yFssS)
   % R = x(1); % Resistance in Ohms
   % L = x(2); % Inductance in Henrys
    L = x(1); % Inductance in Henrys
    L1 = x(2); % Inductance in Henrys

    % Impedance of RLC circuit (example)
    %Z = abs((((((1./(C*((w1)^2)))+((1./(C*((w1)^2)))/(((w2)^2/(w1)^2-1))))./(C)-w.^2*(1./(C*((w1)^2)))*((1./(C*((w1)^2)))/(((w2)^2/(w1)^2-1))))./(1j*w*(1./(C*((w1)^2)))+1./(1j*w*C))).^(-1)));
    yFssC = abs((((L+L1)./(C)-w.^2*L*L1)./(1j*w*L+1./(1j*w*C))).^(-1));
    %Z(isinf(Z)) = 40000;
    cost = sum(abs(yFssC - yFssS).^2)
    % figure
    % plot(w, yFssC)
    % hold on
    % plot(w, yFssS)
    % hold off
end