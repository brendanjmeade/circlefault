close all;

nu = 0.25;
nseg = 1;

npts = 50;
xvec = linspace(-2, 2, npts);
yvec = linspace(-2, 2, npts);
[xmat, ymat] = meshgrid(xvec, yvec);
xobs = xmat(:);
yobs = ymat(:);

dip = 90;
ldep = 1;
bdep = 0;


% Create a bunch of elements in circle
[cx1, cy1, cx2, cy2] = discretizedarc(0, 360, 1, 100);
uxtotal = zeros(numel(xobs), 1);
uytotal = zeros(numel(xobs), 1);
uztotal = zeros(numel(xobs), 1);

figure; hold on;
for i=1:numel(cx1)
    [strike, L, W, ofx, ofy, ~, ~, ~, ~, ~, ~] = fault_params_to_okada_form(cx1(i), cy1(i), cx2(i), cy2(i), deg_to_rad(dip), ldep, bdep);
    [ux, uy, uz] = okada_plus_op(ofx, ofy, strike, ldep, deg_to_rad(dip), L, W, 1, 0, 0, xobs, yobs, nu);
    uxtotal = uxtotal + ux;
    uytotal = uytotal + uy;
    uztotal = uztotal + uz;
end

quiver(xobs, yobs, uxtotal, uytotal)
% contourf(xmat, ymat, reshape(uztotal, npts, npts))

for i=1:numel(cx1)
    plot([cx1(i), cx2(i)], [cy1(i), cy2(i)], "-k")
end
axis equal;
box on;
xlim([-2, 2])
ylim([-2, 2])



function [x1, y1, x2, y2] = discretizedarc(thetastart, thetaend, radius, npts)
% Create geometry of discretized arc
    thetarange = linspace(thetastart, thetaend, npts+1);
    x = radius * cosd(thetarange);
    y = radius * sind(thetarange);
    x1 = x(1:1:end-1);
    x2 = x(2:1:end);
    y1 = y(1:1:end-1);
    y2 = y(2:1:end);
end
