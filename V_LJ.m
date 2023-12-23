function V = V_LJ(r, A, B, r_c)

    V_c = A*(r_c^(-12)) - B*(r_c^(-6));
    V = zeros(size(r));
    V( r <= r_c) = A*(r(r <= r_c).^(-12)) - B*(r(r <= r_c).^(-6)) - V_c;
    
end