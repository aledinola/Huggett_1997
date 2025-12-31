function [V1] = sub_howard(V1,apol_ind,payoff,pi_z,beta,n_howard)

[na,nz] = size(V1);

payoff_max = zeros(na,na);
for z_c = 1:nz
    for a_c = 1:na
        payoff_max(a_c,z_c) = payoff(apol_ind(a_c,z_c),a_c,z_c);
    end
end

for h_c = 1:n_howard % Howard iteration counter
    for z_c = 1:nz
        %for a_c = 1:na
        %EVh = dot(z_prob(z_c,:),V1(kpol_ind(a_c,z_c),:));
        EVh = V1(apol_ind(:,z_c),:)*pi_z(z_c,:)';
        V1(:,z_c) = payoff_max(:,z_c)+beta*EVh;
        %end
    end
end

end %end function 

