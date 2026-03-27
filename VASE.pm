package VASE;
use strict;
use warnings;
use PDL;
use PDL::Fit::LM;
use PDL::NiceSlice;

# Constructor
sub new {
    my ($class, %args) = @_;
    my $self = {
        layers => $args{layers} || 1,
        model  => $args{model}  || sub { die "Model not implemented" },
        data   => undef,
    };
    bless $self, $class;
    return $self;
}

# Load data from file
sub load_data {
    my ($self, $filename) = @_;
    open my $fh, '<', $filename or die "Cannot open $filename: $!";
    my @lines = <$fh>;
    close $fh;
    my @data;
    for my $line (@lines) {
        next if $line =~ /^\s*#/;   # skip comment lines
        next if $line =~ /^\s*$/;    # skip blank lines
        my @fields = split ' ', $line;
        push @data, \@fields;
    }
    my $data = pdl \@data;   # This will be (n x m) where n is rows, m is columns.
    $self->{data} = $data;
    return $data;
}

# Set model function
sub set_model {
    my ($self, $model) = @_;
    $self->{model} = $model;
}

# Fit data to model
sub fit {
    my ($self, $initial_params) = @_;
    my $data = $self->{data};
    my $model = $self->{model};

    # data shape is (nfields, npts); dim0=fields, dim1=data points
    # Extract x as (npts, 2) so $x->(:,0)=wavelength, $x->(:,1)=angle
    my $x_data = $data->(0:1,:)->xchg(0,1);

    # Build y to match model output order: [psi_all, delta_all]
    my $y_data = $data->((2),:)->flat->append($data->((3),:)->flat);
    my $sigma  = ones($y_data->nelem);

    # lmfit sizes $dyda from $x->getdim(0), so pass a dummy x whose
    # first dimension equals the number of y values (2*npts)
    my $x_fit = sequence($y_data->nelem);

    # Wrapper: adapts user model to lmfit interface ($x, $par, $ym, $dyda)
    my $fit_func = sub {
        my ($x, $par, $ym, $dyda) = @_;

        my $y_model = &$model($par, $x_data);
        $ym .= $y_model;

        # Numerical partial derivatives via finite differences
        my $np  = $par->nelem;
        my $eps = 1e-7;
        for my $i (0 .. $np - 1) {
            my $par_h = $par->copy;
            $par_h->slice("($i)") += $eps;
            $dyda->slice(",($i)") .= (&$model($par_h, $x_data) - $y_model) / $eps;
        }
    };

    my ($ym, $finalp, $covar, $iters) = lmfit(
        $x_fit, $y_data, $sigma, $fit_func, $initial_params,
        {Maxiter => 300, Eps => 1e-7}
    );

    return $finalp;
}

1;
