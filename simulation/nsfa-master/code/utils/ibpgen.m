function Z = ibpgen(N, alpha);
    Z = zeros(N, 0);
    K=0;
    for i = 1:N
        entry_probabilities = sum(Z(1:(i-1),:),1) / i;
        Z(i,:) = rand(size(entry_probabilities)) ...
            < entry_probabilities;
        num_new_columns = poissrnd( alpha/i);
        Z(i,K+1:K+num_new_columns)=1;
        K=K+num_new_columns;
    end