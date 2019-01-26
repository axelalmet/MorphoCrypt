function ContactCondition = HasRodHitRegionContact(contactSol, parameters)
% Function to determine whether or not the rod should transition from point
% contact to region contact, as determined by a (close to) zero bending 
% moment at the contact point. 

x = contactSol.y(2,:);
m = contactSol.y(7,:);
w0 = parameters.w0;
contactIndex = find(x == w0, 1);

ContactCondition = 0;
if (m(contactIndex) < 0)
    
    ContactCondition = 1;
    
end