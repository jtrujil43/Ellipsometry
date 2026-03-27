# Copilot Instructions for Ellipsometry Analysis Codebase

## Overview

This is a Perl-based ellipsometry data analysis toolkit for fitting optical models to experimental ellipsometry measurements. The primary components are a `VASE` module for model fitting using the Levenberg-Marquardt algorithm and analysis scripts for processing spectroscopic ellipsometry data.

## Dependencies & Environment

- **Perl**: Core language
- **PDL (Perl Data Language)**: Scientific computing, matrix operations
  - `PDL::Fit::LM`: Levenberg-Marquardt fitting
  - `PDL::NiceSlice`: Matrix slicing syntax (enables `$x->(:,0)` notation)

## Running Tests

```bash
# Run the main fitting example
perl -I. test_fit.pl

# Note: Must include current directory in @INC with -I. since VASE.pm is local
```

## Architecture

### Core Components

1. **VASE.pm**: Main analysis module
   - Constructor accepts `layers` parameter for multi-layer models
   - `load_data()`: Reads whitespace-delimited data files with comment support
   - `set_model()`: Accepts user-defined fitting functions
   - `fit()`: Performs Levenberg-Marquardt fitting using PDL::Fit::LM

2. **test_fit.pl**: Example implementation demonstrating linear model fitting

### Data Format

Input files follow this column structure:
```
# Wavelength(nm) Angle(deg) Psi(deg) Delta(deg)
400 70 45.0 120.0
410 70 44.5 121.0
```

- Lines starting with `#` are comments and skipped
- Blank lines are ignored
- Data expected as whitespace-separated columns

### Model Function Convention

Model functions receive two parameters:
- `$params`: PDL piddle containing fit parameters
- `$x`: PDL piddle with input data (wavelength in column 0, angle in column 1)

Must return flattened PDL containing concatenated Psi and Delta values:
```perl
sub model {
    my ($params, $x) = @_;
    my $wavelength = $x->(:,0);  # PDL NiceSlice syntax
    
    # Calculate Psi and Delta
    my $psi = ...;
    my $delta = ...;
    
    return cat($psi, $delta)->flat;  # Concatenate and flatten
}
```

## Key Conventions

### PDL NiceSlice Syntax
- `$x->(:,0)` extracts first column (all rows)
- `$data->(:,2:3)->flat` extracts columns 2-3 and flattens
- Requires `use PDL::NiceSlice;` in the module

### Data Access Patterns
- Column 0: Wavelength (nm)
- Column 1: Angle of incidence (degrees)  
- Column 2: Psi (ellipsometry parameter, degrees)
- Column 3: Delta (ellipsometry parameter, degrees)

### Fitting Setup
- Initial parameters passed as PDL piddle: `pdl [a, b, c, d]`
- `lmfit` expects input data in specific format:
  - `x`: Independent variables (wavelength, angle)
  - `y`: Dependent variables (flattened Psi/Delta measurements)
  - `f`: Model function closure

## Data Organization

The `data/` directory contains experimental datasets organized by:
- Equipment: `Jovan_Ellipsometer/`, `OTS_on_SiO2/`
- Material systems: `tantalum oxide/`, `Aluminum Oxide/`
- Sample types and dates in subdirectories

## Development Notes

- Module path resolution requires `perl -I.` for local imports
- Error handling uses die() for file operations and missing models
- All fitting uses Levenberg-Marquardt via PDL::Fit::LM
- Data structures are PDL piddles throughout for numerical efficiency