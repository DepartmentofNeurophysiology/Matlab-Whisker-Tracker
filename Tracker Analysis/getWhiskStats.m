function [Manual,Tracker] = getWhiskStats(m_theta,m_peaks,m_troghs,t_theta,t_peaks,t_troghs)
conflag = zeros(1,length(m_peaks)-1);
proflag = zeros(1,length(m_peaks)-1);
dt = (1/300)*1000;
Tracker.protraction_amplitude = [];
Tracker.contraction_amplitude = [];
Tracker.protraction_duration = [];
Tracker.contraction_duration = [];
Tracker.contraction_speed = [];
Tracker.protraction_speed = [];

Manual.protraction_amplitude = [];
Manual.contraction_amplitude = [];
Manual.protraction_duration = [];
Manual.contraction_duration = [];
Manual.contraction_speed = [];
Manual.protraction_speed = [];

peakidxP = zeros(1, length(m_peaks)-1);
peakidxC = zeros(1, length(m_peaks)-1);

for i = 1:length(m_peaks)-1
    
    
    % Find manual parameters
    Manual.cycle_length(i) = m_peaks(i+1)-m_peaks(i);
    tid = find(m_troghs-m_peaks(i) > 0 ,1,'first');
    if isempty(tid)
        continue
    end
    Manual.protraction_amplitude(i) = m_theta(m_peaks(i)) - m_theta(m_troghs(tid));
    Manual.protraction_duration(i) = dt*(m_troghs(tid)-m_peaks(i));
    Manual.protraction_speed(i) = Manual.protraction_amplitude(i)/(Manual.protraction_duration(i)/1000);
    Manual.contraction_amplitude(i) = m_theta(m_troghs(tid)) - m_theta(m_peaks(i+1));
    Manual.contraction_duration(i) = dt*(m_peaks(i+1)-m_troghs(tid));
    Manual.contraction_speed(i) = Manual.contraction_amplitude(i)/(Manual.contraction_duration(i)/1000);
    
    [~, idx] = min(abs(t_peaks-m_peaks(i)));
    if ~isempty(idx)
        peakidxP(i) = t_peaks(idx);
        tid = find(t_troghs-t_peaks(idx) > 0,1,'first');

        if ~isempty(tid)
            peakidxC(i) = t_troghs(tid);
            Tracker.protraction_amplitude(i) = t_theta(t_peaks(idx)) - t_theta(t_troghs(tid));
            Tracker.protraction_duration(i) = dt*(t_troghs(tid)-t_peaks(idx));
            Tracker.protraction_speed(i) = Tracker.protraction_amplitude(i)/(Tracker.protraction_duration(i)/1000);


            proflag(i) = 1;
            if tid < length(t_peaks) & idx < length(t_peaks)
                Tracker.contraction_duration(i) = dt*(t_peaks(idx+1) - t_troghs(tid));
                Tracker.contraction_amplitude(i) =  t_theta(t_troghs(tid))-t_theta(t_peaks(idx+1));
                Tracker.contraction_speed(i) = Tracker.contraction_amplitude(i)/(Tracker.contraction_duration(i)/1000);
                conflag(i) = 1;
            end
        end
    end
    
    
end

Manual.pro_id = m_peaks(proflag==1);
if length(m_troghs) <= length(conflag)
    Manual.con_id = m_troghs(conflag(1:length(m_troghs))==1);
else
    Manual.con_id = m_troghs(conflag == 1);
end

Tracker.pro_id = peakidxP;
Tracker.con_id = peakidxC;

Manual.protraction_amplitude = Manual.protraction_amplitude(proflag == 1);
Manual.protraction_duration = Manual.protraction_duration(proflag == 1);
Manual.protraction_speed = Manual.protraction_speed(proflag == 1);

Manual.contraction_amplitude = Manual.contraction_amplitude(conflag == 1);
Manual.contraction_duration = Manual.contraction_duration(conflag == 1);
Manual.contraction_speed = Manual.contraction_speed(conflag == 1);

Tracker.protraction_amplitude = Tracker.protraction_amplitude(proflag == 1);
Tracker.protraction_duration = Tracker.protraction_duration(proflag == 1);
Tracker.protraction_speed = Tracker.protraction_speed(proflag == 1);

Tracker.contraction_amplitude = Tracker.contraction_amplitude(conflag == 1);
Tracker.contraction_duration = Tracker.contraction_duration(conflag == 1);
Tracker.contraction_speed = Tracker.contraction_speed(conflag == 1);




end