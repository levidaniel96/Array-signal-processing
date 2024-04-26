function [B] = ULA(psi,w,n_m_n0)

v_theta=exp(1i*(n_m_n0)'*psi);
B=w'*v_theta;
end


