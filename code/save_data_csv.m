% Author: Alastair Garner
% email: alastairgarner@outlook.com
% August 2019; Last revision: 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function save_data_csv(printArray,figureName)

    printArray = printArray';
    fid = fopen(strcat(figureName,'.csv'),'w');
    fprintf(fid,'%s,%s,%s,%s\n',string({'genotype','group','id','value'}));
    fprintf(fid,'%s,%d,%d,%f\n',printArray{:});
    fclose(fid);

end