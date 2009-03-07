function val = get(p, propName)
% GET Get f_hat properties from the specified object
% and return the value
switch propName
case 'Data'
    val = p.f_hat;
otherwise
    error([propName,' Is not a valid f_hat property']);
end
