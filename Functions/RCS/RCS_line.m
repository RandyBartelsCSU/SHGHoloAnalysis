function [D_R,TH] = RCS_line(E_tot,x,z,k,R_test,dR)

E_y = squeeze(E_tot(end/2,:,end));
z_half = z(end);
TH = atan(x./z_half);
R = sqrt(x.^2 + z_half^2);
E_R = E_y;%.*(R);
D_R = abs(E_R).^2;
D_R = D_R./max(D_R);

end

