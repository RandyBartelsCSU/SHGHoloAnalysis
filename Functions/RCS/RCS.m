function [D_R,TH] = RCS(E_tot,X,Z,k,R_test,dR)
% Given a total field, calculate the RCS of the scattered field

E_yz = squeeze(E_tot(end/2,:,:));
Z_renorm = Z-Z(1,1);
E_inc = E_yz(1,1)*exp(1i*k*Z_renorm);
E_sca = E_yz-E_inc.';
[TH,R,E_R] = cart2pol(X,Z,E_sca);
TH(R<(R_test-dR)| R>(R_test+dR))=[];
E_R(R<(R_test-dR)| R>(R_test+dR))=[];
R(R<(R_test-dR)| R>(R_test+dR))=[];
E_R(TH<0 | TH>pi/2)=[];
R(TH<0 | TH>pi/2)=[];
TH(TH<0 | TH>pi/2)=[];
[TH,I_TH] = sort(TH);
E_R=E_R(I_TH);
R=R(I_TH);
D_R = abs(E_R).^2;
D_R = D_R./max(D_R);

end

