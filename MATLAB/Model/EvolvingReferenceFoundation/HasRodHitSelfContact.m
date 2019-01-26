function SelfContact = HasRodHitSelfContact(rodSol, parameters)
% Function to determine whether or not rod has hit self-contact, which is
% determined by the condition x = x0, where x0 is a given parameter.

w0 = parameters.w0;
theta = rodSol.y(6,:);
x = rodSol.y(2,:);
contactIndices = find(abs(theta + 0.5*pi) < 1e-3);

SelfContact = 0;

if (min(x(contactIndices)) < w0)
    
    SelfContact = 1;
    
end
