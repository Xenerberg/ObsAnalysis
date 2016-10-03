clc;clear all;
%Define rotational components
syms q_0 q_1 q_2 q_3 om_x om_y om_z ita_1 ita_2 ita_3 ita_0 p_x p_y p_z del_q_1 del_q_2 del_q_3 del_q_0 del_ita_1 del_ita_2 del_ita_3 del_ita_0;
%Define translational components
syms r_x r_y r_z r_dot_x r_dot_y r_dot_z rho_t_x rho_t_y rho_t_z;
q_v = [q_1;q_2;q_3];
q = [q_v;q_0];
ita_v = [ita_1;ita_2;ita_3];
ita = [ita_v;ita_0];
omega = [om_x;om_y;om_z];
p = [p_x;p_y;p_z];
r = [r_x; r_y; r_z];
r_dot = [r_dot_x;r_dot_y;r_dot_z];
rho_t = [rho_t_x;rho_t_y;rho_t_z];

%Define state
X = [q;omega;p;r;r_dot;rho_t;ita];
%X = [q;omega;p;r;r_dot];
%Define vector fields
psi_p_w = [p(1)*om_y*om_z;p(2)*om_x*om_z;p(3)*om_x*om_y];
q_dot = 0.5*fn_CrossTensor([omega;0],0)*q;
omega_dot = psi_p_w;
p_dot = zeros(3,1);
r_dot = r_dot;
r_ddot = zeros(3,1);
rho_t_dot = zeros(3,1);
ita_dot = zeros(4,1);
f = [q_dot;omega_dot;p_dot;r_dot;r_ddot;rho_t_dot;ita_dot];
%f = [q_dot;omega_dot;p_dot;r_dot;r_ddot];


%Define observation equation
h_1 = r + fn_CreateRotationMatrix(q)*rho_t; %rho_c = 0 for now
h_2 = fn_CrossTensor(ita,0)*q;
h_3 = q.'*q - 1;
h_4 = ita.'*ita - 1;
h = [h_1;h_2;h_3;h_4];
%h = [h_1;h_2];
row_1 = jacobian(h,X);

L_1 = simplify(row_1*f);
row_2 = jacobian(L_1,X);

L_2 = simplify(row_2*f);
row_3 = jacobian(L_2,X);

L_3 = simplify(row_3*f);
row_4 = jacobian(L_3,X);

L_4 = simplify(row_4*f);
row_5 = jacobian(L_4,X);

O = [row_1;row_2;row_3;row_4];

O_val = double(subs(O,{q_0,q_1,q_2,q_3,om_x,om_y,om_z,p_x,p_y,p_z,r_x,r_y,r_z,r_dot_x,r_dot_y,r_dot_z,rho_t_x,rho_t_y,rho_t_z,ita_0,ita_1,ita_2,ita_3},{a(4),a(1),a(2),a(3) .1,.1,.1, .1,.1,.1, .1,.2,.1, .1,.1,.1, 0.1,0.1,0.1, 1,0,0,0}));

rank(double(subs(row_1,{q_0,q_1,q_2,q_3,om_x,om_y,om_z,p_x,p_y,p_z,r_x,r_y,r_z,r_dot_x,r_dot_y,r_dot_z,rho_t_x,rho_t_y,rho_t_z,ita_0,ita_1,ita_2,ita_3},{1,0,0,0, .1,.1,.1, .1,.1,.1, .1,.1,.1, .1,.1,.1, 0.0,0.0,0.0, 1,0,0,0})))
	

%fn_CreatePsi(q,rho_t);
jacobian(L_1(1:3),q)

psi = fn_CreatePsi(q,rho_t);
psi_q = jacobian(psi*q_dot,q);
psi_om = psi*fn_CrossTensor(q,1);
psi_rho_t = jacobian(psi*q_dot,rho_t);

theta_1 = psi_q*q_dot;
theta_2 = psi*fn_CrossTensor(q,1)*[psi_p_w;0]

L_2_h1 = theta_1 + theta_2;

theta_q = jacobian(L_2_h1,q);
theta_om = jacobian(L_2_h1,omega);
theta_p = jacobian(L_2_h1,p);
theta_rho_t = jacobian(L_2_h1,rho_t);

