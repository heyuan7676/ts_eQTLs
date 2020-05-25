function Y=sample_data(param_set, N, D)

Y = param_set.G * param_set.X + mvnrnd(zeros(N,D),diag(param_set.lambdae.^-1))';