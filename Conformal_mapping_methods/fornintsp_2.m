%%% FORNBERGS METHOD %%%

function [f,s,erri] = fornintsp_2(n,itmax,x,y)
% Fornberg's method for interior regions from B. Fornberg, A numerical
% method for conformal mappings, SIAM J. Sci. Stat. Comput., 1 (1980)
% 386-400 or T. K. DeLillo and J. A. Pfaltzgraff, Numerical conformal
% mapping methods for simply and doubly connected regions, SIAM J. Sci.
% Comput., 19 (1998) 155-171.
%
% n = number of Fourier points.
% itmax = number of iterations.
% tl = length of curve parmeter interval.
% curve = a cell array containing the information about a boundary
% curve. The first element is the name of a boundary curve. The rest
% are parameters which determine the boundary curve. If the curve is
% a spline, the last element of the cell is a case for the spline.
% guess = a cell array determining the initial guess for boundary ... correspondence.
% center = f(0).
% f(n) = f(s(n)) n complex points along curve at s.
% s(n) = n values of parameter 0 ≤ s(n) ≤ tl.
% erri = successive iteration error.
% errfx = discretization error for conformal mapping.
% errsx = discretization error for boundary correspondence.
% See also WEGAI, WEGSM, WDRHCG.s = tl*(k - 1)/n;
% Thomas K. DeLillo, Lianju Wang 06-29-2000.
n2 = n/2; center=0;
[x1, x2, x3, y1, y2, y3, h, tl]  = spline_(x,y); % fit x,y with spline
s = tl*(0:n-1)'/n; % Initial guess for boundary correspondence.
% Start of Fornberg iteration loop.
for it = 1:itmax
    [f,e] = interp_1(x,x1,x2,x3,y,y1,y2,y3,h,tl,s);
    ne = abs(e); e = e./ne; c = fft(f);
    % Apply I(-,N) = diag(1,0,...,0,1,...,1).
    c(2:n2+1) = 0;
    c = ifft(c);
    a1 = real(center*conj(e)); % Fix f(0) = center.
    b = -real(c.*conj(e)) + a1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% CONJUGATE GRADIENT METHOD %%%%%%
    % initialization.
    v = zeros(n,1); r = b; p = r; rr =dot(r,r);
    % The start of conjugate gradient iteration.
    for in = 1:n
        if (norm(p) < 1.0e-14) break, end
        c = e.*p;
        c = fft(c);
        % Apply I(-,N) = diag(1,0,...,0,1,...,1).
        c(2:n2+1) = 0; c = ifft(c); q = real(c.*conj(e));
        q(1) = q(1) + p(1)*n/4; % Add Q*p, Q rank 1 to fix f(1).
        alpha = rr/dot(p,q); v = v + alpha*p;
        r = r - alpha*q; rr1 = dot(r,r);
        if (sqrt(rr1) < 1.0e-14) break,end
        beta = rr1/rr; p = r + beta*p; rr = rr1;
    end
    %%%%%% end of conjugate gradient iteration %%%%%%
    sl = s; s = s + v./ne; erri(it) = norm(s-sl,inf);
end
% end of outer iteration.
erri = erri(:);
end

%%% PARAMETERIZING THE BOUNDARY %%%
function [x1,x2,x3,y1,y2,y3,h,tl] = spline_(x,y)
x = x(:); y = y(:);
if abs(x(1) - x(end)) > 100*eps | abs(y(1) - y(end)) > 100*eps
  x(end+1) = x(1);
  y(end+1) = y(1);
end
n1 = length(x); n = n1 - 1;
dx = diff(x); dy = diff(y);
h = sqrt(dx.^2 + dy.^2); tl = sum(h); h(n1) = h(1);
p = h(1:n); q = h(2:n1);
a = q./(p + q); b = 1 - a;
c = spdiags([[b(n);ones(n-1,1)] [a(2:n);0] [2*ones(n,1)] [0;b(1:n-1)] ...
    [ones(n-1,1);a(1)]] ...
            , [-n+1 -1 0 1 n-1], n, n);
d1 = 3*(a.*dx./p + b.*[dx(2:n); x(2) - x(n1)]./q);
mmdflag = spparms('autommd');
spparms('autommd', 0);
x1 = c\d1;
spparms('autommd', mmdflag);
x1(2:n1) = x1; x1(1) = x1(n1);
d = 3*(a.*dy./p + b.*[dy(2:n); y(2) - y(n1)]./q);
mmdflag = spparms('autommd');
spparms('autommd', 0);
y1 = c\d;
spparms('autommd', mmdflag);
y1(2:n1) = y1; y1(1) = y1(n1);
x2(2:n1) = 2*(x1(1:n) + 2*x1(2:n1) - 3*dx./p)./p;
y2(2:n1) = 2*(y1(1:n) + 2*y1(2:n1) - 3*dy./p)./p;
x2(1) = x2(n1); y2(1) = y2(n1);
x2 = x2'; y2 = y2';
x3 = diff(x2)./p; x3(n1) = x3(1); y3 = diff(y2)./p; y3(n1) = y3(1);
end


%%% INTERPOLATION %%%
function [fval,fder] = interp_1(x,x1,x2,x3,y,y1,y2,y3,h,tl,s)
x = x(:); y = y(:); s = s(:);
n1 = length(x); nn = n1 - 1; n = length(s);
s = s - floor(s/tl)*tl;
cumsumh = [0; cumsum(h(1:nn))];
ppx = mkpp(cumsumh, [x3(1:nn)/6; x2(1:nn)/2; x1(1:nn); x(1:nn)]);
ppy = mkpp(cumsumh, [y3(1:nn)/6; y2(1:nn)/2; y1(1:nn); y(1:nn)]);
fval = ppval(ppx,s) + i*ppval(ppy,s);
ppdx = mkpp(cumsumh, [x3(1:nn)/2; x2(1:nn); x1(1:nn)]);
ppdy = mkpp(cumsumh, [y3(1:nn)/2; y2(1:nn); y1(1:nn)]);
fder = ppval(ppdx,s) + i*ppval(ppdy,s);
end

% projection method - Wegmann survey p. 390
N = 512; alpha=.8; iter = 4000; t = 2*pi*[0:N-1]/N; S = t;
eta = @(S) cos(S) + i*alpha*sin(S);
eta_dot = @(S) -sin(S) + i*alpha*cos(S);
for k=1:iter
  etas=eta(S);
  B = fft(etas)/N;
  Bn = [2*B(1) i*imag(B(2)) zeros(1,N/2-1) i*imag(B(N/2+2)) 2*B(N/2+3:N)];
  gk = ifft(Bn);
  etad = eta_dot(S);
  U = - real(gk./etad);
  norm(U,inf)
  S = S + U;
end
  etas = eta(S);
  B = fft(etas)/N;
  Bp = [0 real(B(2)) B(3:N/2+1) real(B(N/2+2)) zeros(1,N/2-2)];
  etat = eta(t);
  plot(etas,'.')
  hold on;
  plot(etat,'x');
  axis equal
% Error =max(abs(etat-etas))

% Plot the real and imaginary parts of f(s)
figure;
subplot(2,1,1);
plot(s, real(f), 'b', s, imag(f), 'r');
xlabel('Parameter s');
ylabel('Real and Imaginary Parts');
title('Real and Imaginary Parts of f(s)');
legend('Real(f(s))', 'Imag(f(s))');
grid on;

% Plot the magnitude and phase of f(s)
% subplot(2,1,2);
% plot(s, abs(f), 'g', s, angle(f), 'm');
% xlabel('Parameter s');
% ylabel('Magnitude and Phase');
% title('Magnitude and Phase of f(s)');
% legend('|f(s)|', 'Phase(f(s))');
% grid on;

% Plot the boundary curve
figure;
plot(x, y, 'k');
xlabel('x');
ylabel('y');
title('Boundary Curve');
axis equal;
grid on;


% Define parameters of the ellipse
a = 2; % Semi-major axis length
b = 1; % Semi-minor axis length
theta = linspace(0, 2*pi, 100); % Angles for parameterization

% Parametric equations for an ellipse centered at the origin in the complex plane
z = a * exp(1i * theta); % Complex representation of the ellipse

% Assign the real and imaginary parts to x and y
x = real(z);
y = imag(z);

% Call the fornintsp2 function with the ellipse coordinates
n = 1000; % Number of Fourier points
itmax = 10; % Number of iterations
[f, s, erri] = fornintsp_2(n, itmax, x, y);

disp(f);
