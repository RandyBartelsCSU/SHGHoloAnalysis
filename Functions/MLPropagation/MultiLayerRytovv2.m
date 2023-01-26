function evol=MultiLayerRytovv2(fxx,fyy,lambda,n_imm,dz,V,U_in,Eps,opt)


% This function computes the propagation of an initial field U_in through a
% scattering potential V using a first Rytov approximation for each slice of V. 

% Inputs: fxx,fyy: Spatial frequency grids (ifftshifted)
%         lambda: wavelength in microns
%         n_imm: index of refraction of background material
%         dz:    Thickness of layer in microns
%         V:     Scattering potential of object/material 3d matrix array
%         U_in:  input field
%         eps:   absorption term
%         opt:   options either 'out' or 'Vol'. 'out' just outputs the
%         final output field. 'Vol' outs a 3d array containing the field as it propagates through each
%         slice

SquareRt=@(a) abs(real(sqrt(a)))+1i*abs(imag(sqrt(a))); % making sure imaginary and real part is pos otherwise will get gain or backward propagating field

prop_phs= 1i*2*pi*SquareRt((n_imm/lambda)^2-(fxx.^2+fyy.^2)); %Add small absorption term to avoid indeterminates in the angular greens function

prop=@(z) exp(prop_phs*z);
%Mask=(n_imm/lambda)^2>1.01*((fxx.^2+fyy.^2));

AG= ((-1i.*exp(prop_phs.*dz)./(4.*pi.*SquareRt((n_imm/lambda)^2-(fxx.^2+fyy.^2)+1i*Eps)))); % Angular Greens function
%AG(isnan(AG)==1)=0;
propdz=prop(dz);
U=U_in; % Initial Field
switch opt
    case 'out'
%S=zeros(size(V));
for i=1:size(V,3)
    S=U;
    U=ifft2(propdz.*(fft2(U))); % Prop dz to next layer
    
    Us=ifft2((fft2(S.*V(:,:,i))).*(AG))*dz; % Scattered field of nth layer
    %Us2=ifft2((fft2(Us.*V(:,:,i))).*(AngGreensn(Mask,fxx,fyy,lambda,n_imm,prop_phs,dz)))*dz;
    %U=U+Us;
    %U=U.*exp((Us./U)+(Us2./U)-.5*(Us2./U).^2);
    B1=Us./U;
    B1(abs(B1)>1)=0; % Avoiding indeterminants/large numbers (from zeros in field) assumes phase is small
%     B2=Us2./U;
%     B2(abs(B2)>1)=0; % Added these to avoid indeterminants/large numbers
    %U=U.*exp((B1)+(B2)-.5*((B2)).^2);
    U=U.*exp((B1));
end
evol=U;

    case 'Vol'
        S=zeros(size(V));
for i=1:size(V,3)
    S(:,:,i)=U;
    U=ifft2(propdz.*(fft2(U))); % Prop dz to next layer
    
    Us=ifft2((fft2(S(:,:,i).*V(:,:,i))).*(AG))*dz; % Scattered field of nth layer
    %Us2=ifft2((fft2(Us.*V(:,:,i))).*(AngGreensn(Mask,fxx,fyy,lambda,n_imm,prop_phs,dz)))*dz;
    %U=U+Us;
    %U=U.*exp((Us./U)+(Us2./U)-.5*(Us2./U).^2);
    B1=Us./U;
    B1(abs(B1)>1)=0;% Avoiding indeterminants/large numbers (from zeros in field) assumes phase is small
%     B2=Us2./U;
%     B2(abs(B2)>1)=0; % Added these to avoid indeterminants/large numbers
    %U=U.*exp((B1)+(B2)-.5*((B2)).^2);
    U=U.*exp((B1));
end
S(:,:,end+1)=U;
evol=S;

end

evol=evol(:,:,2:end);


% function AG=AngGreensn(Mask,fxx,fyy,lambda,n_imm,prop_phs,z)
% 
% AG= ((-1i.*exp(prop_phs.*z)./(4.*pi.*sqrt((n_imm/lambda)^2-(fxx.^2+fyy.^2))))).*Mask;
% AG(isnan(AG)==1)=0;
% end
end