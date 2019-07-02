function f = LinearInterpolation(f0, f1, timeSteps)
theta1 = (0:timeSteps) / timeSteps;
theta0 = (timeSteps:-1:0) / timeSteps; 
f = f1' * theta1 + f0' * theta0;
end