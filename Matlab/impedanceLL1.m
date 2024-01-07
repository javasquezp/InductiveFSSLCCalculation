% Equivalent Circuit Model Function
function Y = impedanceC(x, C, w)
   % R = x(1); % Resistance in Ohms
   % L = x(2); % Inductance in Henrys
    L = x(1); % Inductance in Henrys
    L1 = x(2); % Inductance in Henrys

    % Impedance of RLC circuit (example)
    %Z = abs((((((1./(C*((w1)^2)))+((1./(C*((w1)^2)))/(((w2)^2/(w1)^2-1))))./(C)-w.^2*(1./(C*((w1)^2)))*((1./(C*((w1)^2)))/(((w2)^2/(w1)^2-1))))./(1j*w*(1./(C*((w1)^2)))+1./(1j*w*C))).^(-1)));
    Y = abs(((1e3*(L+L1)./(C)-1e3*w.^2*L*L1)./(1j*w*L+1./(1j*w*C))).^(-1));
    %Z(isinf(Z)) = 40000;

end