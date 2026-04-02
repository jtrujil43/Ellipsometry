use strict;
use warnings;
use PDL;
use PDL::NiceSlice;
use Physics::Ellipsometry::VASE;
use PDL::Constants qw(PI E);

# Create VASE object with 1 layer
my $vase = Physics::Ellipsometry::VASE->new(layers => 1);

# Load sample data
$vase->load_data('w1_11012006.dat');

# Define model function (linear model example)
sub model {
    my ($params, $x) = @_;
    
    # Unpack parameters
    my $a = $params->(0);
    my $b = $params->(1);
    my $c = $params->(2);
    my $d = $params->(3);
    
    # Compute linear model (using only wavelength)
    my $wavelength = $x->(:,0);   # first column: wavelength
    # Note: $x->(:,1) contains angle (available for more complex models)

    my $psi = $a - $b * $wavelength;
    my $delta = $c + $d * $wavelength;
    
    return cat($psi, $delta)->flat;
}

sub cauchy_model {
	my ($params, $x) = @_;
	
	# Unpack parameters
	my $a = $params->(0);
	my $b = $params->(1);
	my $c = $params->(2);
	
	my $n0 = $params->(3); # ambient index
	my $n2 = $params->(4); # substrate index
	my $d = $params->(5); # thickness [nm]
	
	# 
	# Unpack independent vars
	#
	my $lambda = $x->(:,0); # wavelength [nm]
	my $theta0 = $x->(:,1) * (PI / 180.0); # incident angle [radians]
	
	#
	# Film refractive index from Cauchy
	#
	my $n1 = $a + $b / ($lambda**2) + $c / ($lambda**4);
	my $k1 = zeroes($lambda); # transparent
	
	#
	# Snell's law
	#
	my $sin_theta1 = $n0 * sin($theta0) / $n1;
	my $theta1 = asin($sin_theta1);
	
	my $sin_theta2 = $n0 * sin($theta0) / $n2;
	my $theta2 = asin($sin_theta2);
	
	#
	# Phase thickness beta
	# beta = (2*PI / lambda) n1 d cos theata1
	#
	my $beta = (2 * PI / $lambda) * $n1 * $d * cos($theta1);
	
	# 
	# Fresnel coefficients
	#
	
	# Air/film
	my $r01s = ($n0*cos($theta0) - $n1*cos($theta1))
             / ($n0*cos($theta0) + $n1*cos($theta1));
             
	my $r01p = ($n1*cos($theta0) - $n0*cos($theta1))
		     / ($n1*cos($theta0) + $n0*cos($theta1));
		     
	# Film/substrate
	my $r12s = ($n1*cos($theta1) - $n2*cos($theta2))
			 / ($n1*cos($theta1) + $n2*cos($theta2));
			 
	my $r12p = ($n2*cos($theta1) - $n1*cos($theta2))
	         / ($n2*cos($theta1) + $n1*cos($theta2));
	         
	# 
	# Thin-film Fresnel reflectances
	# r = (r01 + r12 exp(-2 i beta)) / (1 + r01 r12 exp(-2 i beta))
	#

	my $phase = exp(-2*i*$beta);
	my $rs = ($r01s + $r12s*$phase) / (1 + $r01s*$r12s*$phase);
	my $rp = ($r01p + $r12p*$phase) / (1 + $r01p*$r12p*$phase);

	#
	# Ellipsometric ratio
	#
	my $rho = $rp / $rs;

	# 
	# Psi and Delta
	#

	my $psi = atan( abs($rho) ); # tan(psi) = |rp/rs|
	my $delta = arg($rho);		 # delata = phase(rp/rs)

	return cat($psi, $delta)->flat;
	
}

# $vase->set_model(\&model);
$vase->set_model(\&cauchy_model);

# Initial parameters: [a, b, c, d] for linear model (exact solution for sample data)
# my $initial_params = pdl [65, 0.05, 80, 0.1];
my $initial_params = pdl [65, 0.05, 80, 1.0, 4.0, 1.0];

# Perform fit
my $fit_params = $vase->fit($initial_params);

# Extract results
my ($a, $b, $c, $n0, $n2, $d) = list $fit_params;
print "Fit results:\n";
print "a: $a\n";
print "b: $b\n";
print "c: $c\n";
print "n0: $n0\n";
print "n2: $n2\n";
print "d: $d\n";
print "  Psi   = $a - $b * wavelength\n";
print "  Delta = $c + $d * wavelength\n";
