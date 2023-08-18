function [n_events, line_scan_avg_proc, event_int, event_loc] = find_Ca_events(segmented_cells, wave_thresh, temporal_res)
%{ 
Identify calcium (Ca) events from line scan images of individual cells.

INPUTS: 

-segmented_cells: cell array of 2D numerical matrices representing portions
of the original line scan image containing cells. Each cell array element
contains the image for an individual cell. The number of elements in the
cell array will be equal to the number of cells found. 

-wave_thresh: percent of max mean signal needed to be considered a Ca event

-temporal_res: seconds/pixel, temporal resolution of line scan images


OUTPUTS:

-n_events: a 1D numerical array containing the number of Ca events found in
each image contained in 'segmented_cells'.

-line_scan_avg_proc: an abbreviation of "line scan average processed". This
variable is a 1D cell array where each element contains a 1D numerical
array of the processed time-series data from the cell's line scan. This is
the signal used to identify Ca events who can be visualized by plotting
this variable along side the 'event_int' and 'event_loc' variables.

-event_int: cell array whose elements are 1D numerical arrays of intensity
values for each Ca event identified. The cell array will have a number of
elements equal to 'n_cells' and each element will have a length equal to
the number of Ca events found for that cell. These intensity values are
sampled from the same element in the 'line_scan_avg_proc' variable.

-event_loc: cell array whose elements are 1D numerical arrays of locations
(or time-of-occurance) for each Ca event identified. The cell array will
have a number of elements equal to 'n_cells' and each element will have a
length equal to the number of Ca events found for that cell. These
event locations are sampled from the same element in the
'line_scan_avg_proc' variable.


EXAMPLE: 

If you want to make a line plot of the 1D time-series of Ca intensity for
cell 'i' with dots marking the Ca events, you could use the following code:

i = 1; % change as needed to whatever cell you want to plot for
hold on; % tell Matlab you want to render multiple plots on the same graph
plot(line_scan_avg_proc{i},'k'); % line graph of 1D Ca intensity (y) over time (x)
scatter(event_loc{i_cell}, event_int{i_cell},'filled','r'); % plot a dot for each Ca event
%}

    n_cells = numel(segmented_cells);

    % preallocate output variables
    n_events = zeros(1,n_cells);
    line_scan_avg_proc = cell(1,n_cells);
    event_int = cell(1,n_cells);
    event_loc = cell(1,n_cells);

    for i_cell = 1:n_cells

        % process temporal signal to identify Ca events
        cell_line_scan = segmented_cells{i_cell};
        cell_line_scan = im2double(cell_line_scan); % convert to double so we can do math
        line_scan_avg_raw = sum(cell_line_scan, 2); % condense scan to 1D
        line_scan_avg_proc{i_cell} = line_scan_avg_raw - mean(cell_line_scan, 'all'); % subtract off mean to diminish noise
        line_scan_avg_proc{i_cell}(line_scan_avg_proc{i_cell} < 0) = 0; % delete negative signal
        window_width_sec = 0.1024; % seconds
        window_width_pix = round(window_width_sec / temporal_res); % number of pixels for moving average
        line_scan_avg_proc{i_cell} = movmean(line_scan_avg_proc{i_cell}, window_width_pix); % smooth signal using moving mean with 50-pixel window
        line_scan_avg_proc{i_cell} = line_scan_avg_proc{i_cell} - min(line_scan_avg_proc{i_cell}); % subtract min
        line_scan_avg_proc{i_cell} = line_scan_avg_proc{i_cell} / max(line_scan_avg_proc{i_cell}); % normalize to [0,1]
    
        % Identify Ca waves and sparks
        [event_int{i_cell}, event_loc{i_cell}] = findpeaks(line_scan_avg_proc{i_cell}, 'minpeakheight', wave_thresh,'MinPeakProminence',wave_thresh); % find all events
        n_events(i_cell) = numel(event_int{i_cell}); % count waves
    end

end