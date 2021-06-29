fucntion [reshaped_strains] = read_exper_CSV(file_name_csv)

x = csvread(file_name_csv);                                                 % read data from CSV file                       
                                
reshaped_strains = mySort(x);                                               % reshape columns with coordinates into organization according to element number

function sorted = mySort(x) 
        xSorted = sortrows(x, 1);
        ySorted = sortrows(xSorted, 2);
        incrementing = 1:size(ySorted(:,3:5), 1);
        sorted = [ transpose(incrementing) ySorted(:,3:5) ] ;    
end  