% Postprocessing
fprintf('Nonlinear static test for Element %d \n',Element)
fprintf(' n & DOFs & ux & uy & uz \n')
for k=1:size(Results,1) 
  fprintf('%d & %d & %10.8f & %10.8f & %10.8f  \n',Results(k,1:5))
end
