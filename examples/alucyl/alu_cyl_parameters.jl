function alu_cyl_parameters()
	# Material parameters of the solid cylinder
	
	E = 70000*phun("MPa")::FFlt;
	nu = 0.33::FFlt;
	rho = 2700*phun("KG/M^3")::FFlt;
	radius = 0.5*phun("ft"); 
	leng = 2*phun("ft"); 
	mass_shift = (2*pi*100) ^ 2; # to resolve rigid body modes
    
	return E, nu, rho, radius, leng, mass_shift
end
