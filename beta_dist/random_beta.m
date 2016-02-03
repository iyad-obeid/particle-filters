function Y = random_beta(mu, var_x)

alpha = (mu.^2 - mu.^3 - mu.*var_x)./var_x;
beta = (mu - 2*mu.^2 + mu.^3 - var_x + mu.*var_x)./var_x;

Y = betarnd(alpha, beta);

