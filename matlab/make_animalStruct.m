%% make_animalStruct

function newStruct = make_animalStruct(oldStruct)

fields2change = {'animal','time','crabspeed','speed','x','y','unique_numb'};

% oldStruct = dStruct;

for ii = 1:length(oldStruct)
    nums = [oldStruct(ii).unique_numb];
%     unique(nums)
    lengths = accumarray(nums',1);
    lengths = nonzeros(lengths);
    
    for jj = 1:length(fields2change)
        val = oldStruct(ii).(fields2change{jj});
        oldStruct(ii).(fields2change{jj}) = mat2cell(val,1,lengths);
    end
end

newStruct = oldStruct;
