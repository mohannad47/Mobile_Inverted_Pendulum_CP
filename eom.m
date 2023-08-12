function qdd = eom(params, th, phi, dth, dphi, u)
  % This is the starter file for the week5 assignment

  % Provided params are
  % params.g: gravitational constant
  % params.mr: mass of the "rod"
  % params.ir: rotational inertia of the rod
  % params.d: distance of rod CoM from the wheel axis
  % params.r: wheel radius

  % Provided states are:
  % th: wheel angle (relative to body)
  % phi: body pitch
  % dth, dphi: time-derivatives of above
  % u: torque applied at the wheel
g = params.g  ;
m = params.mr ;
i = params.ir  ;
l = params.d;
r = params.r ;


A = -m*r^2;
B = (-m*r^2 - m*r*l*cos(phi));
C = -m*r*l*(sin(phi)*dphi^2)-u;

A2 = m*r^2+m*r*l*cos(phi);
B2= -(-m * r^2 - 2*m*r*l*cos(phi) - m*l^2-i);
C2 = -(-m*r*l*sin(phi)*dphi^2-m*g*l*sin(phi));
  
AA = [ A B ; A2 B2];
BB = [C; C2]; 
XX = linsolve(AA , BB);
thdd = XX(1);
phidd = XX(2);

qdd = [thdd;phidd] ;
  % THE STUDENT WILL FILL THIS OUT
end