% Ca event test script. Analyze one cell at a time.

close all;
clear all;
clc;

%% Adjustable parameters

cell_thresh = 0.2; % percent of total average signal needed to be considered a cell
wave_thresh = 0.1; % percent of max mean signal needed to be considered a Ca2+ wave
temporal_res = 1.9531/1000;
spatial_res = 0.21; % microns/pixel, spatial resolution of line scan images

%% Import Image

[file, path] = uigetfile('*.tif');
line_scan_image = imread(fullfile(path,file));

%% Find Ca Events

[segmented_cells, n_cells, cell_edges, avg_int_timeAx] = find_cells(line_scan_image, cell_thresh, spatial_res);
[n_events, line_scan_avg_proc, event_int, event_loc] = find_Ca_events(segmented_cells, wave_thresh, temporal_res);

%% Plot Results

[n_time,n_space] = size(line_scan_image);

% plot origianl line scan with indivudal cells outlined
figure;
subplot(1,2,1);
hold on;
plot(-avg_int_timeAx,0:(numel(avg_int_timeAx)-1),'k');
plot([-cell_thresh,-cell_thresh],[0,numel(avg_int_timeAx)]);
ylim([0,numel(avg_int_timeAx)]);
legend(["Time Average Ca Trace","Cell Threshold"]);
title("Time Average");
set(gca,'xticklabel',[],'yticklabel',[]);
hold off;
subplot(1,2,2);
hold on;
imagesc(line_scan_image'); % plot full line scan
for i_cell = 1:n_cells
    plot([0,n_time],[cell_edges{i_cell}(1),cell_edges{i_cell}(1)],'color','r','linewidth',1); % plot cell upper outline
    p = plot([0,n_time],[cell_edges{i_cell}(2),cell_edges{i_cell}(2)],'color','r','linewidth',1); % plot cell lower outline
end
legend(p,"Cell Outlines");
title(sprintf("Full Line Scan (%i Cells)",n_cells));
xlim([0,n_time]);
padding = 10;
ymax = n_space + (2 * padding);
ylim([-padding,ymax]);
hold off;
set(gca,'xticklabel',[],'yticklabel',[]);

% plot individual cell line scans, Ca trace, and Ca events
for i_cell = 1:n_cells
    figure; % tell Matlab to make a new figure, one for each cell
    sgtitle("Cell " + num2str(i_cell));
    subplot(2,1,1);
    imagesc(segmented_cells{i_cell}');
    title("Line Scan");
    set(gca,'xticklabel',[],'yticklabel',[]);
    subplot(2,1,2);
    hold on; % tell Matlab you want to render multiple plots on the same graph
    plot(line_scan_avg_proc{i_cell},'k'); % line graph of 1D Ca intensity (y) over time (x)
    scatter(event_loc{i_cell}, event_int{i_cell},'filled','r'); % plot a dot for each Ca event
    title("Ca^{2+} Trace " + "(" + num2str(n_events(i_cell)) + " events)");
    xlabel("Time");
    ylabel("Ca Intensity");
    legend(["Ca Trace","Ca Events"]);
    hold off;
end
