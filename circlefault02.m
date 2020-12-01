close all;
clear all;

% plot_schematic()
% basic_profiles()
% rotelasticinter()
% rdgrids()
% approx_scaling()

function approx_scaling()
    nu = 0.25;
    nseg = 100;
    dip = 90;
    D = 25;
    currentr = 200;
    rcurrentr = 100;
    npts = 200;
    xvec = linspace(0.001, rcurrentr, npts);
    yvec = linspace(0, rcurrentr, npts);
    [xmat, ymat] = meshgrid(xvec, yvec);
    xobs = xmat(:);
    yobs = ymat(:);

    % Create circle
    [rcx1, rcy1, rcx2, rcy2] = discretizedarc(0, 360, rcurrentr, 100);
    [cx1, cy1, cx2, cy2] = discretizedarc(0, 360, currentr, 100);
    cx1 = cx1 + currentr;
    cx2 = cx2 + currentr;
    uxshallow = zeros(numel(xobs), 1);
    uyshallow = zeros(numel(yobs), 1);
    uxrot = zeros(numel(xobs), 1);
    uyrot = zeros(numel(yobs), 1);

    % rotation component
    for i=1:length(xobs)
        if inpolygon(xobs(i), yobs(i), cx1, cy1)
            temp = cross([xobs(i)-currentr, yobs(i), 0], [0, 0, 1/currentr]);
            uxrot(i) = temp(1);
            uyrot(i) = temp(2);
        end
    end

    % Elastic component
    for i=1:numel(cx1)
        [strike, L, W, ofx, ofy, ~, ~, ~, ~, ~, ~] = fault_params_to_okada_form(cx1(i), cy1(i), ...
                                                          cx2(i), cy2(i), ...
                                                          deg_to_rad(dip), 20, 0);
        [ux, uy, ~] = okada_plus_op(ofx, ofy, strike, 20, deg_to_rad(dip), L, W, 1, 0, 0, xobs, yobs, nu);
        uxshallow = uxshallow + ux;
        uyshallow = uyshallow + uy;
    end
    umagshallow = sqrt(uxshallow.^2 + uyshallow.^2);
    umagrot = sqrt(uxrot.^2 + uyrot.^2);
    umagtotal = sqrt((uxrot - uxshallow).^2 + (uyrot - uyshallow).^2);


    % Small finite source
    [strike, L, W, ofx, ofy, ~, ~, ~, ~, ~, ~] = fault_params_to_okada_form(0, -5, ...
                                                      0, 5, ...
                                                      deg_to_rad(dip), 20, 0);
    [uxf, uyf, ~] = okada_plus_op(ofx, ofy, strike, 20, deg_to_rad(dip), L, W, 1, 0, 0, xobs, yobs, nu);
    umagf = sqrt(uxf.^2 + uyf.^2);
    umagf = 1-umagf;

    % Calculate displacment as a function of distance from the 
    dist = sqrt(xobs.^2 + yobs.^2);
    for i=1:length(xobs)
        if ~inpolygon(xobs(i), yobs(i), rcx1, rcy1)
            umagtotal(i) = 0;
            umagf(i) = 0;
            dist(i) = 0;
        end
    end
    delidx = find(dist == 0);
    dist(delidx) = [];
    umagtotal(delidx) = [];
    umagf(delidx) = [];

    % Remove observations that are not in circle
    [dist, distsortidx] = sort(dist);
    umagtotal = umagtotal(distsortidx);
    umagf = umagf(distsortidx);

    figure("Position", [0, 0, 600, 300]);
    hold on;
    plot(linspace(0, rcurrentr, 1000), 0.5+1/pi*atan(linspace(0, rcurrentr, 1000)/D), '-k', "linewidth", 2)
    sp = csaps(dist, umagf, 0.5);
    fnplt(sp, '-b');
    sp = csaps(dist, umagtotal, 0.0005);
    fnplt(sp, '-r');
    lh = legend("infinitely long", "finite", "circle");
    set(lh, "fontsize", 12);
    set(lh, "location", "southeast");
    legend boxoff;
    box on;
    xlabel("distance (km)");
    ylabel("v (mm/yr)");
    xlim([0.0, 50.0])
    xticks([0.0, 25, 50])
    ylim([0.25, 1.0])
    yticks([0.25, 0.5, 0.75, 1.0])
    yticklabels({"0.25", "0.50", "0.75", "1.00"})
    set(gca, "TickDir", "out");
    axis square;
    set(gcf, "color", "w");
    export_fig("approx_scaling.pdf");
end

function rdgrids()
% Loop over radius and locking depth
    nu = 0.25;
    nseg = 100;
    r = 150; % radius of circle
    dip = 90;
    ldep = 1e14;
    bdep = 20;
    nprofpts = 10000;

    ngrid = 20;
    xprof = linspace(-200, 200, nprofpts);
    rvec = linspace(10, 100, ngrid);
    Dvec = linspace(0, 25, ngrid);
    [rmat, Dmat] = meshgrid(rvec, Dvec);
    rmatflat = rmat(:);
    Dmatflat = Dmat(:);
    umaxvalmat = zeros(size(rmatflat));
    umaxlocmat = zeros(size(rmatflat));
    for i=1:length(rmatflat)
        i / length(rmatflat) * 100
        uytotal = circle_fault_profile(xprof, rmatflat(i), Dmatflat(i));
        [maxval, maxidx] = max(uytotal);
        umaxvalmat(i) = maxval;
        umaxlocmat(i) = xprof(maxidx) / rmatflat(i);
    end

    figure("Position", [0, 0, 600, 300]);
    set(gcf, "color", "w")
    subplot(1, 2, 1);
    contourf(rmat, Dmat, reshape(umaxvalmat, ngrid, ngrid), 10);
    ch = colorbar;
    ch.Ticks = [0, 0.5, 1.0]; %Create 8 ticks from zero to 1
    ch.TickLabels = {"0.0", "0.5", "1.0"};
    ch.YLabel.String = "v/v_0";
    caxis([0 1]);
    axis square;
    xlabel("r (km)");
    ylabel("d (km)");
    set(gca,'TickDir','out');
    box on;

    subplot(1, 2, 2);
    contourf(rmat, Dmat, reshape(umaxlocmat, ngrid, ngrid), 10);
    caxis([0.0, 0.25]);
    ch = colorbar;
    ch.Ticks = [0, 0.05, 0.10, 0.15, 0.20, 0.25]; %Create 8 ticks from zero to 1
    ch.TickLabels = {"0.0", "0.05", "0.10", "0.15", "0.20", "0.25"};
    ch.YLabel.String = "x_{max} / r";
    axis square;
    xlabel("r (km)");
    ylabel("d (km)");
    set(gca,'TickDir','out');
    box on;
end



function basic_profiles()
% Loop over radius
    nu = 0.25;
    nseg = 100;
    r = 150; % radius of circle
    dip = 90;
    ldep = 2000000;
    bdep = 20;
    nprofpts = 1000;
    xprof = linspace(-200, 200, nprofpts);

    rvec = [25, 50, 100, 250, 500];
    figure('Position', [0 0 300 300])
    set(gcf, "color", "w")
    hold on;
    plot(xprof, 1/pi*atan((xprof)./bdep)+0.5, ":k", "linewidth", 2);
    for i=1:length(rvec)
        uytotal = circle_fault_profile(xprof, rvec(i), bdep);
        plot(xprof, uytotal, "linewidth", 2);
    end
    legend("SB73",...
           sprintf("r = %d (km)", rvec(1)), sprintf("r = %d (km)", rvec(2)),...
           sprintf("r = %d (km)", rvec(3)), sprintf("r = %d (km)", rvec(4)),...
           sprintf("r = %d (km)", rvec(5)), "location", "northwest")
    yticks([-0.5, 0.0, 0.5 1.0]);
    yticklabels({"-0.5", "0.0", "0.5", "1.0"});
    xlabel("x (km)");
    ylabel("v (mm/yr)");
    set(gca,'TickDir','out');
    box on;
    print("-dpng", "-r300", "basic_profiles.png");
    export_fig("basic_profiles.pdf");
end



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


function uytotal = circle_fault_profile(xprof, r, D)
% Create a bunch of elements in circle
    nseg = 100;
    nprofpts = 1000;
    ldep = 1e10;
    dip = 90;
    nu = 0.25;
    
    [cx1, cy1, cx2, cy2] = discretizedarc(0, 360, r, nseg);
    cx1 = cx1 + r;
    cx2 = cx2 + r;

    % xprof = linspace(-200, 200, nprofpts);
    yprof = zeros(size(xprof));

    % Start of profiles
    uytotal = zeros(numel(xprof), 1);
    for i=1:numel(cx1)
        [strike, L, W, ofx, ofy, ~, ~, ~, ~, ~, ~] = fault_params_to_okada_form(cx1(i), cy1(i), cx2(i), cy2(i), deg_to_rad(dip), ldep, D);
        [~, uy, ~] = okada_plus_op(ofx, ofy, strike, ldep, deg_to_rad(dip), L, W, 1, 0, 0, xprof, yprof, nu);
        uytotal = uytotal + uy;
    end
end


function plot_schematic()
    nseg = 100;
    r = 150; % radius of circle
    nprofpts = 1000;
    xmin = -300;
    xmax = 300;
    ymin = -300;
    ymax = 300;

    xprof = linspace(-300, 300, nprofpts);
    yprof = zeros(size(xprof));

    % Create a bunch of elements in circle
    [cx1, cy1, cx2, cy2] = discretizedarc(0, 360, r, nseg);
    cx1 = cx1 + r;
    cx2 = cx2 + r;

    npts = 20;
    xvec = linspace(0.9*xmin, 0.9*xmax, npts);
    yvec = linspace(0.9*ymin, 0.9*ymax, npts);
    [xmat, ymat] = meshgrid(xvec, yvec);
    xobs = xmat(:);
    yobs = ymat(:);

    % Cartoon block velocities
    uxtotal2d = zeros(numel(xobs), 1);
    uytotal2d = zeros(numel(xobs), 1);
    uxtotalcirc = zeros(numel(xobs), 1);
    uytotalcirc = zeros(numel(xobs), 1);
    for i=1:numel(xobs)
        if xobs(i) > 0
            uytotal2d(i) = 1;
        end
        if inpolygon(xobs(i), yobs(i), cx1, cy1)
            temp = cross([xobs(i)-r, yobs(i), 0], [0, 0, 1/r]);
            uxtotalcirc(i) = temp(1);
            uytotalcirc(i) = temp(2);
        end
    end

    figure('Position', [0 0 600 300])
    set(gcf, "color", "w")
    vel_scale = 20;
    subplot(1, 2, 1);
    hold on;
    fh = fill([0, xmax, xmax, 0], [ymin, ymin, ymax, ymax], "y");
    set(fh, "facecolor", 0.9*[1, 1, 1]);
    plot(xprof, zeros(size(xprof)), "-r", 'linewidth', 2);
    plot([0, 0], [xmin, xmax], "-k", "linewidth", 3.0);
    quiver(xobs, yobs, vel_scale*uxtotal2d, vel_scale*uytotal2d, 0, 'color', [0 0 1]);
    xlabel("x (km)");
    ylabel("y (km)");
    axis equal;
    xlim([xmin, xmax])
    ylim([xmin, xmax])
    set(gca,'TickDir','out');
    box on;

    subplot(1, 2, 2);
    hold on;
    fh = fill(cx1, cy1, "y");
    set(fh, "facecolor", 0.9*[1, 1, 1]);
    plot(xprof, zeros(size(xprof)), '-r', 'linewidth', 2);
    for i=1:numel(cx1)
        plot([cx1(i), cx2(i)], [cy1(i), cy2(i)], "-k", 'linewidth', 3)
    end
    quiver(xobs, yobs, vel_scale*uxtotalcirc, vel_scale*uytotalcirc, 0, 'color', [0 0 1]);
    xlabel("x (km)");
    ylabel("y (km)");
    axis equal;
    xlim([xmin, xmax])
    ylim([xmin, xmax])
    set(gca,'TickDir','out');
    box on;
    print("-dpng", "-r300", "schematic.png");
    export_fig("schematic.pdf");
end


function rotelasticinter()
% Experiments with small blocks
    nu = 0.25;
    nseg = 100;
    dip = 90;
    D = 15;
    currentr = 50;
    npts = 100;
    xvec = linspace(-100, 100, npts);
    yvec = linspace(-100, 100, npts);
    [xmat, ymat] = meshgrid(xvec, yvec);
    xobs = xmat(:);
    yobs = ymat(:);

    % Create circle
    [cx1, cy1, cx2, cy2] = discretizedarc(0, 360, currentr, 100);
    cx1 = cx1 + currentr;
    cx2 = cx2 + currentr;
    uxshallow = zeros(numel(xobs), 1);
    uyshallow = zeros(numel(yobs), 1);
    uxrot = zeros(numel(xobs), 1);
    uyrot = zeros(numel(yobs), 1);

    % rotation component
    for i=1:length(xobs)
        if inpolygon(xobs(i), yobs(i), cx1, cy1)
            temp = cross([xobs(i)-currentr, yobs(i), 0], [0, 0, 1/currentr]);
            uxrot(i) = temp(1);
            uyrot(i) = temp(2);
        end
    end

    % Elastic component
    for i=1:numel(cx1)
        [strike, L, W, ofx, ofy, ~, ~, ~, ~, ~, ~] = fault_params_to_okada_form(cx1(i), cy1(i), ...
                                                          cx2(i), cy2(i), ...
                                                          deg_to_rad(dip), 20, 0);
        [ux, uy, ~] = okada_plus_op(ofx, ofy, strike, 20, deg_to_rad(dip), L, W, 1, 0, 0, xobs, yobs, nu);
        uxshallow = uxshallow + ux;
        uyshallow = uyshallow + uy;
    end
    umagshallow = sqrt(uxshallow.^2 + uyshallow.^2);
    umagrot = sqrt(uxrot.^2 + uyrot.^2);
    umagtotal = sqrt((uxrot - uxshallow).^2 + (uyrot - uyshallow).^2);

    figure("Position", [0, 0, 900, 300]);
    % contourvec = 0:0.1:1;
    contourvec = linspace(0, 1, 100);
    subplot(1, 3, 1)
    hold on;
    contourf(xmat, ymat, reshape(umagrot, npts, npts), contourvec, 'edgecolor','none')
    for i=1:numel(cx1)
        plot([cx1(i), cx2(i)], [cy1(i), cy2(i)], "-k", "linewidth", 1)
    end
    axis equal;
    box on;
    xticks([-100, 0, 100]);
    yticks([-100, 0, 100]);
    xlabel("x (km)");
    ylabel("y (km)");
    ch = colorbar;
    colormap(cool);
    ch.Ticks = [0, 0.5, 1.0]; %Create 8 ticks from zero to 1
    ch.TickLabels = {"0.0", "0.5", "1.0"};
    caxis([0 1]);
    title("rotation", "fontweight", "normal")

    subplot(1, 3, 2)
    hold on;
    contourf(xmat, ymat, reshape(umagshallow, npts, npts), contourvec, 'edgecolor','none')
    for i=1:numel(cx1)
        plot([cx1(i), cx2(i)], [cy1(i), cy2(i)], "-k", "linewidth", 1)
    end
    axis equal;
    box on;
    xticks([-100, 0, 100]);
    yticks([-100, 0, 100]);
    xlabel("x (km)");
    ylabel("y (km)");
    ch = colorbar;
    colormap(cool);
    ch.Ticks = [0, 0.5, 1.0]; %Create 8 ticks from zero to 1
    ch.TickLabels = {"0.0", "0.5", "1.0"};
    caxis([0 1]);
    title("elastic", "fontweight", "normal")

    subplot(1, 3, 3)
    hold on;
    contourf(xmat, ymat, reshape(umagtotal, npts, npts), contourvec, 'edgecolor','none')
    for i=1:numel(cx1)
        plot([cx1(i), cx2(i)], [cy1(i), cy2(i)], "-k", "linewidth", 1)
    end
    axis equal;
    box on;
    xticks([-100, 0, 100]);
    yticks([-100, 0, 100]);
    xlabel("x (km)");
    ylabel("y (km)");
    ch = colorbar;
    colormap(cool);
    ch.Ticks = [0, 0.5, 1.0]; %Create 8 ticks from zero to 1
    ch.TickLabels = {"0.0", "0.5", "1.0"};
    caxis([0 1]);
    title("interseismic", "fontweight", "normal")
    
    export_fig("rotelasticinter.pdf")
end