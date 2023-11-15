function rms = rmsdiff(T1,T2,center)

% implement Jenkinson RMS motion estimate

R = 80; % radius (mm)
M = T2 / T1 - eye(4);
A = M(1:3,1:3);
t = M(1:3,4);
%Tr = trace(A); % this performs sum(trace())
%Tr = diag(A);

rms = sqrt(1/5*R^2 * trace(A' * A) + (t + A*center)' * (t + A*center));

return