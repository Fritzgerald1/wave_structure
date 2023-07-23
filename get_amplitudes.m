function [Amp ] = get_amplitudes(root,wd,Fun)
Amp = zeros(4,length(root));

for ii = 1:length(root)
    if root(1,ii) ~= 0
	   [~,~,A] = Fun(root(ii),wd);
       Amp(:,ii) = A;
    else
       break
    end
end