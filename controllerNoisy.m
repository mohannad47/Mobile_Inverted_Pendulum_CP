
function u = controllerNoisy(params, t, obs)
  % Now you only receive noisy measurements for phi, and must use your EKF  to filter the data and get an estimate of the state
  % obs = [ay; az; gx] with a* in units of g's, and gx in units of rad/s

  % This template code calls the function EKFupdate that you must complete below
  xhat = EKFupdate(params, t, obs);
  phi = xhat(1);
  phidot = xhat(2);

  % The rest of this function should ideally be identical to your solution in week 4
  % Student completes this
persistent t_old
  persistent ei
  if isempty(t_old)
    % initialize
    t_old = 0;
  end
  if isempty(ei)
    % initialize
    ei = 0;
  end
  dt = t - t_old ;
  t_old = t;
  
  ei = ei + dt * phi;
  u = 0.1 *phi + 0.1*phidot+0.5*ei;
end

function xhatOut = EKFupdate(params, t, z)
  % z = [ay; az; gx] with a* in units of g's, and gx in units of rad/s
  %implement a single predict-update step of the EKF
  % Recall ( that you can use persistent variables to create/update any additional state that you need.
  
  persistent t_old
  if isempty(t_old)
    % initialize
    t_old = 0;
  end
  
    persistent x_old
  if isempty(x_old)
    % initialize
    x_old = [0 ; 0];
  end
  
  persistent P_old
  if isempty(P_old)
    % initialize
    P_old =  [1 0 ; 0 10]*0.01;
  end
  
   persistent H
  if isempty(H)
    % initialize
    H = [cosd(x_old(1)) 0; -sind(x_old(1)) 0; 0 1];
  end
  
  persistent h
  if isempty(h)
    % initialize
    h = [sind(x_old(1)) ;cosd(x_old(1)); x_old(2)];

  end

  t_new = t;

      
  A = [1, t_new-t_old ; 0 1];
  Q_tune =  [100 0 ; 0 3]*0.01;
  R_tune = [1 0 0; 0 0.7 0; 0 0 0.01];
  
  x1 = A *  x_old;
  p1 = A * P_old * A' + Q_tune;
  
  K_k = p1 * H' * (H * p1 * H' + R_tune);
  
  xhat = x1 + K_k*((z - h)) ;

  p2 = (eye(2) - K_k * H) * p1;
  
  H = [cosd(xhat(1)) 0; -sind(xhat(1)) 0; 0 1];
  h= [sind(xhat(1)) ;cosd(xhat(1)); xhat(2)];
  
  t_old = t_new;
  x_old = xhat;
  P_old = p2;



  xhatOut = xhat;
end
