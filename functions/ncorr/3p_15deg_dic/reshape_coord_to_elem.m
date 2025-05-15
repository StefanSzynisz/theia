x = csvread("3p_15deg_dic_1.csv")
%x = csvread("3p_15deg_dic_1_missing_1.csv")
%x = csvread("3p_15deg_dic_1_missing_2.csv")
%x = csvread("3p_15deg_dic_1_missing_3.csv")
%x = csvread("3p_15deg_dic_1_missing_4.csv")
%x = csvread("3p_15deg_dic_1_different_1.csv")
%x = csvread("3p_15deg_dic_1_different_2.csv")
%x = csvread("3p_15deg_dic_1_different_3.csv")
%x = csvread("3p_15deg_dic_1_different_4.csv")
%x = csvread("3p_15deg_dic_1_different_5.csv")

function [sorted, nx, ny, delta_x, delta_y] = reshape_coord_to_elem(x)
  % This verifies the input has no missing elements and that elements are
  % evenly spaced
  % It also returns: 
  %  - sorted - the array which has been sorted and with first two columns
  %    removed
  %  - nx, ny - the number of elements on x and y axes
  %  - delta_x, delta_y - the spacing between elements on x and y axes

	x_sorted = sortrows(x, 1)
	y_sorted = sortrows(x_sorted, 2)

  % Need to verify that data is regular grid and that there are no missing elems
  diff_sorted = diff(y_sorted)
  delta_x = diff_sorted(1, 1)
  delta_y_start = diff_sorted(1, 2)
  delta_y = -1
  nx = -1
  ny = -1
  nx_cnt = 0
  ny_cnt = -1
  for i = 1:size(diff_sorted,1)
    if (diff_sorted(i, 1) == delta_x) && (diff_sorted(i, 2) == delta_y_start)
      % Normal elements in the middle of lines
      assert(delta_x == diff_sorted(i, 1)) % Make sure delta_y is the same
      nx_cnt = nx_cnt + 1
    elseif (ny_cnt == -1)
      % Beginning of second y-axis, find out what delta_y is here
      assert((delta_x * nx_cnt) == diff_sorted(i, 1) * -1) % Make sure delta_x is the same
      ny_cnt = 1
      delta_y = diff_sorted(i, 2)
      nx = nx_cnt + 1
      delta_x
      nx
      nx_cnt = 0
    else
      % Normal beginning of new line
      assert((delta_x * nx_cnt) == diff_sorted(i, 1) * -1) % Make sure delta_x is the same
      assert(delta_y == diff_sorted(i, 2)) % Make sure delta_y is the same
      ny_cnt = ny_cnt + 1
      assert((nx_cnt+1) == nx)
      nx_cnt = 0
    end
  end
  assert(nx_cnt == nx - 1) % Make sure we finish on a whole row
  ny = ny_cnt + 1

	incrementing = 1:size(y_sorted(:,3:5), 1)
	sorted = [ transpose(incrementing) y_sorted(:,3:5) ]
end


[sorted, nx, ny, delta_x, delta_y] = reshape_coord_to_elem(x)
