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

    # Levenberg-Marquardt fitting
    my $fit = lmfit {
        p     => $initial_params,
        f     => sub {
            my ($p, $x) = @_;
            my $y_model = $model->($p, $x);
            return $y_model;
        },
        x     => $data->(:,0:1),  # columns 0 and 1: wavelength and angle
        y     => $data->(:,2:3)->flat,  # flattened Psi and Delta
    };

    return $fit->{p};
}

1;
