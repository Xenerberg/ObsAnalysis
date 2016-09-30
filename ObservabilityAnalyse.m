clc;clear all;
syms q_0 q_1 q_2 q_3 om_x om_y om_z ita_1 ita_2 ita_3 ita_0 p_x p_y p_z del_q_1 del_q_2 del_q_3 del_q_0 del_ita_1 del_ita_2 del_ita_3 del_ita_0;
q_v = [q_1;q_2;q_3];
q = [q_v;q_0];
ita_v = [ita_1;ita_2;ita_3];
ita = [ita_v;ita_0];
omega = [om_x;om_y;om_z];
%p = [0.2;0.6;0.1];
p = [p_x;p_y;p_z];
del_q = [del_q_1;del_q_2;del_q_3];
del_ita = [del_ita_1;del_ita_2;del_ita_3];

x = [del_q;omega;p;del_ita];
%x = [del_q;omega;del_ita];
%x = [q;omega;p];

%R = (2*q_0^2 - 1)*eye(3,3) + 2*q_0*fn_VectorToSkewSymmetricTensor(q_v) + 2*q_v*q_v.';
%rho_t = [rho_t_x;rho_t_y;rho_t_z];
psi = [p(1)*om_y*om_z;p(2)*om_x*om_z;p(3)*om_x*om_y];
q_dot = 0.5*fn_CrossTensor([omega;0],0)*q;
omega_dot = psi;
ita_dot = zeros(4,1);
p_dot = zeros(3,1);
del_q_dot = 0.5*fn_VectorToSkewSymmetricTensor(omega)*del_q;
del_ita_dot = zeros(3,1);
f = [del_q_dot;omega_dot;p_dot;del_ita_dot];
%f = [del_q_dot;omega_dot;del_ita_dot];
%f = [q_dot;omega_dot;p_dot];
h = fn_CrossTensor([del_ita;1],0)*[del_q;1];
h = h(1:3);

row_1 = jacobian(h,x);

L_1 = row_1*f;
row_2 = jacobian(L_1,x);

L_2 = row_2*f;
row_3 = jacobian(L_2,x);

L_3 = row_3*f;
row_4 = jacobian(L_3,x);

O = [row_1;row_2;row_3;row_4];


