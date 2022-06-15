function X = tdma(M, D)
%% Classic algo for solving tri-diagonal linear problems
%% Also known as Thomas algo

    N = length(M);
    P = zeros(1, N);
    Q = zeros(1, N);
    X = zeros(1, N);
    A = @(i) M(i,i);
    B = @(i) -M(i, i+1);
    C = @(i) -M(i, i-1);

    P(1) = B(1) / A(1);
    Q(1) = D(1) / A(1);

    for it = 2:N-1
        P(it) = B(it) / (A(it) - C(it) * P(it-1));
        Q(it) = (C(it) * Q(it-1) + D(it)) / (A(it) - C(it)*P(it-1));
    end

    Q(N) = (C(N) * Q(N-1) + D(N)) / (A(N) - C(N)*P(N-1));

    X(end) = Q(end);
    for it = 1:N-1
        itb = N - it;
        X(itb) = P(itb) * X(itb + 1) + Q(itb);
    end
    X = X';
end

