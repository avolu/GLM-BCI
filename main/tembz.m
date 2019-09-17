function [sig] = tembz(sig, param)
%temporally embeds and zscores a signal using the parameters from param struct
aux_emb = sig;
for i=1:param.NumOfEmb
    aux=circshift(aux_sigs, i*param.tau, 1);
    aux(1:2*i,:)=repmat(aux(2*i+1,:),2*i,1);
    aux_emb=[aux_emb aux];
end
% zscore
sig=zscore(aux_emb);
end

