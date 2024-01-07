% Equivalent Circuit Model Function
function Z = impedanceC(x, w1, w2, w)
   % R = x(1); % Resistance in Ohms
   % L = x(2); % Inductance in Henrys
    C = x; % Capacitance in Farads

    % Impedance of RLC circuit (example)
    Z = abs((((((1./(C*((w1)^2)))+((1./(C*((w1)^2)))/(((w2)^2/(w1)^2-1))))./(C)-w.^2*(1./(C*((w1)^2)))*((1./(C*((w1)^2)))/(((w2)^2/(w1)^2-1))))./(1j*w*(1./(C*((w1)^2)))+1./(1j*w*C))).^(-1)));
    %Z(isinf(Z)) = 40000;

end