function ksi = change(x,min_x,max_x);
    fun_ksi = @(x1) (2*x1-(max_x+min_x))/(max_x-min_x);
    ksi=fun_ksi(x);
end