%% format_title

function titlename = format_title(filename)

% filename = fils(1).name;

expr = '\d\d\d\d\d\d\d\d_\d\d\d\d\d\d[@]([\w()]*[@][\w()]*)';

token = regexp(filename,expr,'tokens','once');
names = string(token);

names = strrep(names,'None','');


expr = '([@][a-zA-Z0-9]*)[_](\w*)';
exprrep = '$1^{$2}';
namesfix = regexprep(names,expr,exprrep);

expr = '(^[a-zA-Z0-9]*)([_])(\w*[@])';
exprrep = '$1; $3';
namesfix = regexprep(namesfix,expr,exprrep);

namesfix = strrep(namesfix,'@','>');
namesfix = strrep(namesfix,'_','-');

namesfix = regexprep(namesfix,'SM2>','SM2-null>');
expr = '(SM2)[-]([\w()-]*)([>])';
exprrep = 'dEaat1^{$2}$3';
namesfix = regexprep(namesfix,expr,exprrep);

namesfix = strrep(namesfix,'-NEDA','::Venus');
namesfix = strrep(namesfix,'SM2','dEaat1^{null}');
namesfix = strrep(namesfix,'>S103A','>hEAAT1^{S103A}');
namesfix = strrep(namesfix,'>WT','>hEAAT1^{WT}');
namesfix = strrep(namesfix,'>M128R','>hEAAT1^{M128R}');

titlename = namesfix;

end