# bcalc

Reduced transition probability calculator

Maintainer: Jonathan Williams

## Installation

Use `make` to compile. To run the program from anywhere, move the resulting `bcalc` executable to any directory under your `$PATH` environment variable.

This shouldn't depend on any external libraries.  Tested on CentOS 7.

## Usage

The calculator is run from the command line, for example:

```
bcalc -e VALUE -m VALUE -lt VALUE
```

Running `bcalc` without any parameters will print a list of parameters.

### Parameters

Both of the following are needed:

|**Parameter**|**Description**|
|:---:|:---:|
| -e | transition energy in keV |
| -m |  multipole (eg. E1, M1, E2, etc.) |

One of the following is needed:

|**Parameter**|**Description**|
|:---:|:---:|
| -lt | Mean transition lifetime (in ps) |
| -hl | Transition half-life (in ps) |
| -b | Reduced transition probability (for the L multipole) in units of e^2 fm^(2L) for electric multipoles or uN^2 fm^(2L-2) for magnetic multipoles. |

Optional parameters:

|**Parameter**|**Description**|
|:---:|:---:|
| -d | mixing ratio with the L+1 multipole |
| -ji | inital spin (integer or half-integer) |
| -jf | final spin (integer or half-integer) |

Flags:

|**Flag**|**Description**|
|:---:|:---:|
| --barn | Use/calculate transition probability with spatial dimension in barns rather than fm (eg. B(E2) in e<sup>2</sup> b<sup>2</sup> rather than e<sup>2</sup> fm<sup>4</sup>). |
| --up | Use/calculate transition probability from final to initial state instead of vice versa.  If used, requires the `-ji` and `-jf` parameters. |
| --verbose | Print detailed information. |
| --help | Print a list of parameters. |