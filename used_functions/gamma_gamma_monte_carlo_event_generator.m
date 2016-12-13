function [cos_alpha,cos_theta,theta_3,theta_4,phi_3,phi_4,E_3,E_4,Vq,Ecm_pair,gamma_cm]=gamma_gamma_monte_carlo_event_generator(E_1,Px_1,Py_1,Pz_1,Sx_1,Sy_1,Sz_1,E_2,Px_2,Py_2,Pz_2,Sx_2,Sy_2,Sz_2)
global ev cs cross



%% energy of CM
%Ecm_pair=sqrt(m_1^2+m_2^2+2*E_1*E_2-2*(Px_1*Px_2+Py_1*Py_2+Pz_1*Pz_2));
Ecm_pair=sqrt(2*E_1*E_2-2*(Px_1*Px_2+Py_1*Py_2+Pz_1*Pz_2));

%% cos alpha is interaction angle of two particles in the lab frame
cos_alpha=(Px_1*Px_2+Py_1*Py_2+Pz_1*Pz_2)/(sqrt(Px_1^2+Py_1^2+Pz_1^2)*sqrt(Px_2^2+Py_2^2+Pz_2^2));

%% Generate a theta using monte-carlo rejection_method
% energy of CoM here in MeV;Vq is value of differential crossection
[cos_theta,Vq]=rejection_method(Ecm_pair/1e6);
% theta_out=acos(cos_theta);
%% Generate a phi as uniform on [0 2*pi]
phi=2*pi*rand(1);

%% bild a momentum vector in center of mass system K_3=K_4=Ecm_pair/2;
% K_3x=(Ecm_pair/2)*sin(theta_out)*cos(phi);
% K_3y=(Ecm_pair/2)*sin(theta_out)*sin(phi);
% K_3z=(Ecm_pair/2)*cos(theta_out);
% 
% K_4x=(Ecm_pair/2)*sin(theta_out)*cos(phi+pi);
% K_4y=(Ecm_pair/2)*sin(theta_out)*sin(phi+pi);
% K_4z=(Ecm_pair/2)*cos(pi-theta_out);

K_3x=(Ecm_pair/2)*sqrt(1-cos_theta^2)*cos(phi);
K_3y=(Ecm_pair/2)*sqrt(1-cos_theta^2)*sin(phi);
K_3z=(Ecm_pair/2)*cos_theta;

K_4x=(Ecm_pair/2)*sqrt(1-cos_theta^2)*cos(phi+pi);
K_4y=(Ecm_pair/2)*sqrt(1-cos_theta^2)*sin(phi+pi);
K_4z=-(Ecm_pair/2)*cos_theta;

%% gamma of the center of mass system
gamma_cm = (E_1+E_2)./Ecm_pair;

%% beta_cm
% we put "-" due to Lorenz trasformation mus be P_lab=L(-b)K_cm
bx_cm=-(Px_1+Px_2)/(sqrt(Px_1^2+Py_1^2+Pz_1^2)+sqrt(Px_2^2+Py_2^2+Pz_2^2));
by_cm=-(Py_1+Py_2)/(sqrt(Px_1^2+Py_1^2+Pz_1^2)+sqrt(Px_2^2+Py_2^2+Pz_2^2));
bz_cm=-(Pz_1+Pz_2)/(sqrt(Px_1^2+Py_1^2+Pz_1^2)+sqrt(Px_2^2+Py_2^2+Pz_2^2));

% b_cm=[bx_cm by_cm bz_cm]

b_cm=sqrt(1-1/gamma_cm^2);

%% Lamba matrix of the Lorenz boost
L=[gamma_cm -bx_cm*gamma_cm -by_cm*gamma_cm -bz_cm*gamma_cm; ...
  -bx_cm*gamma_cm 1+(gamma_cm-1)*bx_cm^2/b_cm^2 (gamma_cm-1)*bx_cm*by_cm/b_cm^2 (gamma_cm-1)*bx_cm*bz_cm/b_cm^2; ... 
  -by_cm*gamma_cm (gamma_cm-1)*bx_cm*by_cm/b_cm^2 1+(gamma_cm-1)*by_cm^2/b_cm^2 (gamma_cm-1)*by_cm*bz_cm/b_cm^2; ... 
  -bz_cm*gamma_cm (gamma_cm-1)*bx_cm*bz_cm/b_cm^2 (gamma_cm-1)*by_cm*bz_cm/b_cm^2 1+(gamma_cm-1)*bz_cm^2/b_cm^2;];

%%  transform momentums from CoM to LAB frame 
four_vector_in_lab_3=L*[Ecm_pair/2; K_3x; K_3y; K_3z;];
four_vector_in_lab_4=L*[Ecm_pair/2; K_4x; K_4y; K_4z;];

%% Energy in lab frame
E_3=four_vector_in_lab_3(1);
E_4=four_vector_in_lab_4(1);

%% momentum in lab frame
P_3x=four_vector_in_lab_3(2);
P_3y=four_vector_in_lab_3(3);
P_3z=four_vector_in_lab_3(4);

P_4x=four_vector_in_lab_4(2);
P_4y=four_vector_in_lab_4(3);
P_4z=four_vector_in_lab_4(4);

%% theta in lab frame
theta_3=acos(P_3z/sqrt(P_3x^2+P_3y^2+P_3z^2));
theta_4=acos(P_4z/sqrt(P_4x^2+P_4y^2+P_4z^2));

%% phi in lab frame
phi_3=atan(P_3y/P_3x);
phi_4=atan(P_4y/P_4x);

