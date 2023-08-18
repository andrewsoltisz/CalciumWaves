function [segmented_cells, n_cells, good_edges, avg_int_timeAx] = find_cells(line_scan_image, cell_thresh, spatial_res)
%{ 
Identify and segment individual cells from a line scan image.

INPUTS: 

-line_scan_image: 2D numerical matrix of intesity values composing
the line scan image. 

-cell_thresh: cell segmentation threshold defined as the percent of total
average signal needed to be considered a cell.

-spatial_res: microns/pixel, spatial resolution of the line scan image


OUTPUTS:

-segmeted_cells: cell array of 2D numerical matrices representing portions
of the original line scan image containing cells. Each cell array element
contains the image for an individual cell. The number of elements in the
cell array will be equal to the number of cells found. The cell array will
be empty if no cells are found.

-n_cells: number of cells found in the original line scan image. This
number will be equal to the number of elements composing the
'segmented_cells' output variable.

-cell_edges: cell array of line scan row numbers the segmented cells are
found. Each element in the cell array will be a 2 element numerical array
containing both cell edges. 

-avg_int_timeAx: this is an abbreviation of "average intensity time axis".
This is a 1D numeral array of the intire line scan averaged along the time
axis so as to represent only the spatial component of the data. This
variable is used to identify cells.

%}

    % process spatial signal to identify cells
    avg_int_timeAx = mean(line_scan_image, 1); % average across time domain to see cell signals
    window_width_um = 6.3; % microns
    window_width_pix = round(window_width_um / spatial_res); % number of pixels for moving average 
    avg_int_timeAx = movmean(avg_int_timeAx, window_width_pix); % smooth signal using moving mean
    avg_int_timeAx = avg_int_timeAx - min(avg_int_timeAx); % make sure min signal is = 0
    avg_int_timeAx = avg_int_timeAx ./ max(avg_int_timeAx); % normalize to [0,1]
    avg_int_timeAx = [0,0,avg_int_timeAx,0,0]; % pad with zeros to prevent odd number of derivative peaks

    % find cell edges
    is_a_cell = avg_int_timeAx >= cell_thresh; % identify regions with cells
    cell_edges = find(diff(is_a_cell) ~= 0); % cell edges are marked by sharp changes in intensity
    n_cells = numel(cell_edges) / 2; % each cell as 2 edges

    % create individual images for each cell
    is_whole_num = ((double(int8(n_cells)) - n_cells) == 0);
    segmented_cells = {};
    good_edges = {};
    if n_cells >= 1 && is_whole_num
        i_good_cell = 1;
        for i_cell = 1:n_cells
            idx1 = ((i_cell * 2) - 1);
            idx2 = i_cell * 2;
            crop_width_buffer_um = 2.1; % microns
            crop_width_buffer_pix = round(crop_width_buffer_um / spatial_res); % number of pixels to crop segmented cell image
            cell_edges_1 = cell_edges(idx1) + crop_width_buffer_pix;
            cell_edges_2 = cell_edges(idx2) - crop_width_buffer_pix;
            thresh_width_um = 1.05;
            thresh_width_pix = round(thresh_width_um / spatial_res); % cells must be this wide to pass
            if cell_edges_2 - cell_edges_1 > thresh_width_pix % probably noise/background if signal is this thin
                segmented_cells{i_good_cell} = line_scan_image(:, cell_edges_2:-1:cell_edges_1);
                good_edges{i_good_cell} = [cell_edges_1, cell_edges_2];
                i_good_cell = i_good_cell + 1;
            end
        end
    end
    n_cells = numel(segmented_cells);
end