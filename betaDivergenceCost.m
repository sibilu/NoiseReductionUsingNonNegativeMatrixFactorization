function [bD] = betaDivergenceCost(S,S_ap, beta)
if beta == 0
    bD = sum((S(:)./S_ap(:))-log(S(:)./S_ap(:)) - 1);
elseif beta == 1
    bD = sum(S(:).*(log(S(:))-log(S_ap(:))) + S_ap(:) - S(:));
else
    bD = sum(max(1/(beta*(beta-1))*(S(:).^beta + (beta-1)*S_ap(:).^beta - beta*S(:).*S_ap(:).^(beta-1)),0));
end
end