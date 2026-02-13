function KL1 = KLD(D, A0, A1)

    for i = 1 : length(A0)
        m0(:,i) = A0(i).rnvec(:,1);
        m1(:,i) = A1(i).rnvec(:,1);
        n0(:,i) = A0(i).rnvec(:,2);
        n1(:,i) = A0(i).rnvec(:,2);
    end
    
    m0 = m0';
    m1 = m1';
    n0 = n0';
    n1 = n1';
    mean_m0 = mean(m0);
    mean_m1 = mean(m1);
    
    m0 = cov(m0);
    m1 = cov(m1);
    if det(m0) == 0
        m0 = m0 + 1e-6.*eye(D);
    end
    if det(m1) == 0 
        m1 = m1 + 1e-6.*eye(D);
    end
   
    KL1 = (trace(inv(m1)*m0) + (mean_m1 - mean_m0)*inv(m1)*(mean_m1 - mean_m0)'...
         - D + log(det(m1)/det(m0))) / 2;
     
    mean_n0 = mean(n0);
    mean_n1 = mean(n1);
    
    n0 = cov(n0);
    n1 = cov(n1);
    if det(n0) == 0
        n0 = n0 + 1e-6.*eye(D);
    end
    if det(n1) == 0 
        n1 = n1 + 1e-6.*eye(D);
    end
   
    KL2 = (trace(inv(n1)*n0) + (mean_n1 - mean_n0)*inv(n1)*(mean_n1 - mean_n0)'...
         - D + log(det(n1)/det(n0))) / 2;
     
    KL = (KL1+KL2)/2;
    
end