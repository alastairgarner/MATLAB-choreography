% save_figure.m

% Other m-files required: 
% Subfunctions: 
% MAT-files required: 
% See also:

% Author: Alastair Garner
% email: alastairgarner@outlook.com
% Website: https://github.com/alastairgarner/
% August 2019; Last revision: 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function save_figure(fig_handle,fig_path,fig_type)

    set(fig_handle, 'InvertHardCopy', 'off',...
            'Color', [1 1 1]);

    switch fig_type
        case 'svg'
            print(fig_handle,fig_path,'-dsvg','-painters');
        case 'pdf'
            print(fig_handle,fig_path,'-dpdf','-painters','-fillpage');
        case 'fig'
            savefig(fig_handle,fig_path);
    end
end