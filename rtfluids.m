% Implementation of the seminal papers 'Stable fluids' and 'Real time fluid
% dynamics for games' by Jos Stam. The code includes also alternative
% implementations of the functions 'lin_solve' and 'project' that make use
% of 'gradient' and 'del2': results are in general slightly different due
% to the different way in which edge values are computed.

clear; close all; clc;

%% Demo
n = 128;
dt = 0.1;
diffusivity = 1e-5;
viscosity = 1e-6;
vorticity = 0.3;
force = 5.0;
source = 100.0;

u = zeros(n+2);
v = zeros(n+2);
u0 = zeros(n+2);
v0 = zeros(n+2);
d = zeros(n+2);
d0 = zeros(n+2);

walls = make_walls( ...
    make_box(n, 40, 50, 50, 2) + ...
    make_box(n, 20, 50, 2, 30) + ...
    make_box(n, 20, 78, 50, 2) + ...
    make_box(n, 88, 50, 2, 60));

v0(31,11) = force;
d0(31,11) = source;

map = zeros(size(d));
if ~isempty(walls)
    map(sub2ind(size(map), walls(:,1), walls(:,2))) = 1;
end

figure('Position', [10 10 600 600])
for t = 0:dt:60
    [u, v] = update_velocity(u, v, u0, v0, n, dt, viscosity, vorticity, walls);
    d      = update_density(d, d0, u, v, n, dt, diffusivity, walls);
    show_fluid(u, v, d, 0, map); drawnow;
    if t > 50
        d0 = 0*d0; % stop source
        u0 = 0*u0; % stop force along u
        v0 = 0*v0; % stop force along v
    end
end

function show_fluid(u, v, d, s, map)
    d(map > 0) = 0.5;
    pcolor(d'); colormap(bone); caxis([0 2]); shading(gca, 'interp'); hold on;
    if s; quiver(s*u', s*v', 'AutoScale', 'off', 'Color', 'r'); end
    axis off; hold off
end

function z = make_box(n, x, y, w, h)
    z = zeros(n+2); z(x:x + max(2,w), y:y + max(2,h)) = 1;
end

%% Solver
function x = add_source(x, s, dt)
    x = x + dt*s;
end

function x = set_bounds(x, b, walls)
    if b == 1
        x(  1,2:end-1) = -x(    2,2:end-1);
        x(end,2:end-1) = -x(end-1,2:end-1);
    else
        x(  1,2:end-1) = +x(    2,2:end-1);
        x(end,2:end-1) = +x(end-1,2:end-1);
    end
    if b == 2
        x(2:end-1,  1) = -x(2:end-1,    2);
        x(2:end-1,end) = -x(2:end-1,end-1);
    else
        x(2:end-1,  1) = +x(2:end-1,    2);
        x(2:end-1,end) = +x(2:end-1,end-1);
    end
    x(  1,  1) = 0.5*(x(    2,  1) + x(  1,    2));
    x(  1,end) = 0.5*(x(    2,end) + x(  1,end-1));
    x(end,  1) = 0.5*(x(end-1,  1) + x(end,    2));
    x(end,end) = 0.5*(x(end-1,end) + x(end,end-1));
    if ~isempty(walls)
        x = set_walls(x, b, walls);
    end
end

function w = make_walls(z)
    [w(:,1), w(:,2)] = ind2sub(size(z), find(z > 0));
    z = min(z, 1); ii = 2:length(z)-1; jj = 2:length(z)-1;
    z(ii,jj) = z(ii,jj).*(z(ii-1,jj) + z(ii+1,jj) + z(ii,jj-1) + z(ii,jj+1));
    w(:,3) = z(sub2ind(size(z), w(:,1)-1, w(:,2))) == 0;
    w(:,4) = z(sub2ind(size(z), w(:,1)+1, w(:,2))) == 0;
    w(:,5) = z(sub2ind(size(z), w(:,1), w(:,2)-1)) == 0;
    w(:,6) = z(sub2ind(size(z), w(:,1), w(:,2)+1)) == 0;
end

function x = set_walls(x, b, w)
    if isempty(w); return; end
    q = sum(w(:,3:6), 2); sz = size(x);
    ki = 1; if b == 1; ki = -1; end
    kj = 1; if b == 2; kj = -1; end
    % ki = ones(size(w,1),1); if b == 1; ki(q == 1) = -1; end
    % kj = ones(size(w,1),1); if b == 2; kj(q == 1) = -1; end
    p = + ki.*w(:,3).*x(sub2ind(sz, w(:,1)-1, w(:,2)  )) ...
        + ki.*w(:,4).*x(sub2ind(sz, w(:,1)+1, w(:,2)  )) ...
        + kj.*w(:,5).*x(sub2ind(sz, w(:,1)  , w(:,2)-1)) ...
        + kj.*w(:,6).*x(sub2ind(sz, w(:,1)  , w(:,2)+1));
    ind = sub2ind(sz, w(:,1), w(:,2)); flt = q > 0;
    x(ind(~flt)) = 0; x(ind(flt)) = p(flt)./q(flt);
end

function x = lin_solve(x, x0, n, b, a, c, walls)
    ii = 2:n+1; jj = 2:n+1;
    for iter = 1:20
        x(ii,jj) = (x0(ii,jj) ...
            + a*(x(ii-1,jj) + x(ii+1,jj) + x(ii,jj-1) + x(ii,jj+1)))/c;
        x = set_bounds(x, b, walls);
    end
end

% function x = lin_solve(x, x0, n, b, a, c, walls)
%     ii = 2:n+1; jj = 2:n+1;
%     for iter = 1:20
%         d = 4*del2(x);
%         x(ii,jj) = (x0(ii,jj) + a*(4*x(ii,jj) + d(ii,jj)))/c;
%         x = set_bounds(x, b, walls);
%     end
% end

function x = diffuse(x, x0, n, dt, b, par, walls)
    a = par*dt*n*n; x = lin_solve(x, x0, n, b, a, 1+4*a, walls);
end

function d = advect(d, d0, u, v, n, dt, b, walls)
    dt0 = dt*n; sz = size(d0); ii = 2:n+1; jj = 2:n+1;
    x = ii' - 1 - dt0*u(ii,jj); x = max(0.5, min(x(:), n + 0.5));
    y = jj  - 1 - dt0*v(ii,jj); y = max(0.5, min(y(:), n + 0.5));
    i0 = floor(x) + 1; i1 = i0 + 1; s1 = x - i0 + 1; s0 = 1 - s1;
    j0 = floor(y) + 1; j1 = j0 + 1; t1 = y - j0 + 1; t0 = 1 - t1;
    d(ii,jj) = reshape(...
        + s0.*(t0.*d0(sub2ind(sz,i0,j0))+t1.*d0(sub2ind(sz,i0,j1))) ...
        + s1.*(t0.*d0(sub2ind(sz,i1,j0))+t1.*d0(sub2ind(sz,i1,j1))), n, n);
    d = set_bounds(d, b, walls);
end

function [u, v, p, div] = project(u, v, n, walls)
    ii = 2:n+1; jj = 2:n+1;
    div = zeros(size(u));
    div(ii,jj) = 0.5*n*(u(ii+1,jj)-u(ii-1,jj)+v(ii,jj+1)-v(ii,jj-1));
    div = set_bounds(div, 0, walls);
    p = lin_solve(0*div, -div, n, 0, n*n, 4*n*n, walls);
    u(ii,jj) = u(ii,jj) - 0.5*n*(p(ii+1,jj) - p(ii-1,jj));
    v(ii,jj) = v(ii,jj) - 0.5*n*(p(ii,jj+1) - p(ii,jj-1));
    u = set_bounds(u, 1, walls);
    v = set_bounds(v, 2, walls);
end

% function [u, v, p, div] = project(u, v, n, walls)
%     [~, ux] = gradient(u, 1/n);
%     [vy, ~] = gradient(v, 1/n);
%     div = set_bounds(ux + vy, 0, walls);
%     p = lin_solve(0*div, -div, n, 0, n*n, 4*n*n, walls);
%     [py, px] = gradient(p, 1/n);
%     u = set_bounds(u - px, 1, walls);
%     v = set_bounds(v - py, 2, walls);
% end

function [u, v] = add_vorticity(u, v, n, vorticity, walls)
    ii = 3:n; jj = 3:n;
    fx = zeros(size(u));
    fy = zeros(size(v));
    dvortx = zeros(size(u));
    dvorty = zeros(size(v));
    vort = zeros(size(u));
    vort(ii,jj) = 0.5*n*(v(ii+1,jj)-v(ii-1,jj)-u(ii,jj+1)+u(ii,jj-1));
    dvortx(ii,jj) = 0.5*n*(abs(vort(ii+1,jj))-abs(vort(ii-1,jj)));
    dvorty(ii,jj) = 0.5*n*(abs(vort(ii,jj+1))-abs(vort(ii,jj-1)));
    mag = sqrt(dvortx.^2 + dvorty.^2);
    fx(ii,jj) = +vorticity*dvorty(ii,jj).*vort(ii,jj)./(mag(ii,jj)*n);
    fy(ii,jj) = -vorticity*dvortx(ii,jj).*vort(ii,jj)./(mag(ii,jj)*n);
    fx(mag < 1e-7 | u == 0) = 0;
    fy(mag < 1e-7 | v == 0) = 0;
    fx = set_bounds(fx, 0, walls);
    fy = set_bounds(fy, 0, walls);
    u = u + fx;
    v = v + fy;
    u = set_bounds(u, 1, walls);
    v = set_bounds(v, 2, walls);
end

function x = update_density(x, x0, u, v, n, dt, par, walls)
    x = add_source(x, x0, dt);
    t = x; x = x0; x0 = t; x = diffuse(x, x0, n, dt, 0, par, walls);
    t = x; x = x0; x0 = t; x = advect(x, x0, u, v, n, dt, 0, walls);
end

function [u, v] = update_velocity(u, v, u0, v0, n, dt, par, vorticity, walls)
    u = add_source(u, u0, dt);
    v = add_source(v, v0, dt);
    [u, v] = add_vorticity(u, v, n, vorticity, walls);
    t = u; u = u0; u0 = t; u = diffuse(u, u0, n, dt, 1, par, walls);
    t = v; v = v0; v0 = t; v = diffuse(v, v0, n, dt, 2, par, walls);
    [u, v, u0, v0] = project(u, v, n, walls);
    t = u; u = u0; u0 = t;
    t = v; v = v0; v0 = t;
    u = advect(u, u0, u0, v0, n, dt, 1, walls);
    v = advect(v, v0, u0, v0, n, dt, 2, walls);
    [u, v] = project(u, v, n, walls);
end