
%Set number of simulations to use to critical values
crit_sims=10^6;
alpha_grid=[0.01 0.02 0.05 0.1 0.15 0.2];
gamma_grid=[0.01 0.02 0.05 0.1 0.15 0.2];

%Set up the matrix for output
output=[];

%Sets the matrix for p and k combinations. 
p_grid=[1:1:2];
k_grid=[1:1:50];
p_k_array={p_grid, k_grid};

m=length(p_grid);
p_k_mat=p_k_array(1);
p_k_mat=cell2mat(p_k_mat)';
%%
if m>1
    for n=2:m
        p_k_old=p_k_mat;
        p_k_new=p_k_array(n);
        p_k_new=cell2mat(p_k_new)';
        p_k_mat=[kron(ones(length(p_k_new),1),p_k_old) kron(p_k_new,ones(length(p_k_old),1))];
    end
end

%%
%Calculate value of "a" (the linear combination weight) corresponding to
%gamma_min and different hypotheses
a_min_vec=zeros(length(p_k_mat),1);
crit_vec=zeros(length(p_k_mat),1);
for i=1:length(alpha_grid)
%Set nominal coverage for confidence sets: nominal coverage is 1-alpha
alpha=alpha_grid(i);
    for j=1:length(gamma_grid);
    %Sets the values of gamma_min to be used in two-step confidence set
    %construction
    gamma_min=gamma_grid(j);

        for n=1:length(p_k_mat)
            p=p_k_mat(n,1);
            k=p_k_mat(n,2);

            s = RandStream('mt19937ar','Seed',1);
            RandStream.setGlobalStream(s)
            K_size=sum(randn(crit_sims,p).^2,2);
            J_size=sum(randn(crit_sims,k-p).^2,2);

            crit_quant=@(a) quantile((1+a)*K_size+a*J_size,1-alpha-gamma_min);
            crit_obj=@(a) abs(chi2inv(1-alpha,p)-crit_quant(a));
            a_min_vec(n)=fminsearch(crit_obj,0);

            crit_vec(n)=quantile((1+a_min_vec(n))*K_size+a_min_vec(n)*J_size,1-alpha);
        end
    a = ones(length(p_k_mat),1)*100*alpha;
    g = ones(length(p_k_mat),1)*100*gamma_min;

    output=vertcat(output, horzcat(a,g,p_k_mat,a_min_vec,crit_vec));
    end
end
%%
dlmwrite('a_min.txt',output,'delimiter',' ')
