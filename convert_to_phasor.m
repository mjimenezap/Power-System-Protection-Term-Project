
function voltage_phasor = convert_to_phasor(time_vector, voltage_vector)
V_phasor = zeros(1,size(voltage_vector,1)-80);
    for j=1:(size(voltage_vector,1)-80)

        %Calculating voltage factor
        V_phasor(:,j)=(1/80)*sum(voltage_vector(j+1:j+80,:).*(cos(2*pi*60*time_vector(j+1:j+80,:))+i*sin(2*pi*60*time_vector(j+1:j+80,:)))); 

        coeff_V(1,j)=real(V_phasor(:,j))*sqrt(2); %a1
        coeff_V(2,j)=imag(V_phasor(:,j))*sqrt(2); %a2

        mod_V(:,j) = sqrt(coeff_V(1,j)*coeff_V(1,j) + coeff_V(2,j)*coeff_V(2,j));

        phase_V(:,j) = atan(-coeff_V(2,j)/coeff_V(1,j));

    end
voltage_phasor = [mod_V;phase_V ]
end