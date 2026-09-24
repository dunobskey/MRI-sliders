% 1. Create main window layout
fig = uifigure('Position', [100 100 700 600]); 
ax = uiaxes(fig, 'Position', [60 180 580 390], "LineWidth", 3, 'FontWeight', 'bold');

% 2. Create TR Controls
uilabel(fig, 'Text', 'TR (s):', 'FontSize', 15, 'Position', [40 140 70 22], 'FontWeight', 'bold');
box_tr = uieditfield(fig, 'numeric', 'Position', [115 135 60 25], 'Value', 6.5, 'FontSize', 14, 'FontWeight', 'bold', 'Limits', [0 15]);
sld_tr = uislider(fig, 'Position', [200 150 280 3], 'Limits', [0 15], 'Value', 6.5, 'FontSize', 16, 'FontWeight', 'bold');
sld_tr.MajorTicks = [0,5,10,15];

% 3. Create TE Controls
uilabel(fig, 'Text', 'TE (ms):', 'FontSize', 15, 'Position', [40 85 70 22], 'FontWeight', 'bold');
box_te = uieditfield(fig, 'numeric', 'Position', [115 80 60 25], 'Value', 93, 'FontSize', 14, 'FontWeight', 'bold', 'Limits', [10 400]);
sld_te = uislider(fig, 'Position', [200 95 280 3], 'Limits', [10 400], 'Value', 93, 'FontSize', 16, 'FontWeight', 'bold');

% Package controls into a struct to pass cleanly to callbacks
sliders.tr = sld_tr; 
sliders.te = sld_te; 

% Helper callback functions to keep edit field, slider, and plot synchronized
tr_slider_cb = @(~, val) updateTR(val, box_tr, sld_tr, ax, sliders);
te_slider_cb = @(~, val) updateTE(val, box_te, sld_te, ax, sliders);

% Bind controls
sld_tr.ValueChangedFcn = @(~,~) tr_slider_cb([], sld_tr.Value);
sld_te.ValueChangedFcn = @(~,~) te_slider_cb([], sld_te.Value);

box_tr.ValueChangedFcn = @(~,~) tr_slider_cb([], box_tr.Value);
box_te.ValueChangedFcn = @(~,~) te_slider_cb([], box_te.Value);

% Bind ValueChangingFcn using inline structures to capture continuous live dragging
sld_tr.ValueChangingFcn = @(~,e) tr_slider_cb([], e.Value);
sld_te.ValueChangingFcn = @(~,e) te_slider_cb([], e.Value);

% Initial render
draw(ax, sliders);

function updateTR(val, box, sld, ax, sliders)
    box.Value = val;
    sld.Value = val;
    sliders.tr.Value = val;
    draw(ax, sliders);
end

function updateTE(val, box, sld, ax, sliders)
    box.Value = val;
    sld.Value = val;
    sliders.te.Value = val;
    draw(ax, sliders);
end

function draw(ax, sliders)
    % Extract active timing dimensions from UI controls
    t_r = sliders.tr.Value*1000;
    t_e = sliders.te.Value;
    
    % 1. Tissue Parameters Matrix
    j = [4195, 1800, 1300, 1000; 1885, 150, 100, 80; "CSF", "Blood", "Gray matter", "White matter"];
    R_1relaxivity = 4.71; R_2relaxivity = 5.47;
    C = linspace(0, 1.6, 201)'; % Static contrast baseline range
    
    cla(ax);
    hold(ax, 'on');
    t_i=str2double(j(1,1))*log(2/(1+exp(-t_r/4195)));
    
    % 2. Process math matrix columns across all 4 elements
    for tissue_idx = 1:4
        t_1 = str2double(j(1, tissue_idx));
        t_2 = str2double(j(2, tissue_idx));
        label = j(3, tissue_idx);
        
        M_xy = zeros(201, 1);
        for i = 1:201
            post_t_1 = (1/t_1 + 0.001 * R_1relaxivity * C(i))^-1;
            post_t_2 = (1/t_2 + 0.001 * R_2relaxivity * C(i))^-1;
            x=1-exp(-(t_r-t_i)/post_t_1);
            y=1-(1+x)*exp(-t_i/post_t_1);
            M_xy(i)=y*exp(-t_e/post_t_2);
        end
        plot(ax, C, M_xy, 'LineWidth', 3.5, 'DisplayName', label);
    end
    
    % 3. Format Graphics Window
    grid(ax, 'on');
    legend(ax, 'Location', 'northeast');
    ylim(ax, [0 1]);
    xlim(ax, [0 1.6]);
    xlabel(ax, 'Contrast Agent Concentration (mM)');
    ylabel(ax, 'Signal Intensity M_{xy}/M_0');
    title(ax, sprintf('FLAIR | TR: %.0f ms | TE: %.0f ms | TI: %.0f ms', t_r, t_e, t_i));
    fontsize(ax, 14,'points');
end