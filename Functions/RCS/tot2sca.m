function [E_tot,E_sca] = tot2sca(E_simu,U_inp_end)
%UNTITLED4 Summary of this function goes here

E_tot = squeeze(E_simu(:,:,end));
E_inc = E_simu(1,1,1)*U_inp_end;
E_sca = E_tot-E_inc;

end

