# bcalc

**Reduced transition probability calculator**

Maintainer: Jonathan Williams

## Installation

Use `make` to compile. To run the program from anywhere, move the resulting `bcalc` executable to any directory under your `$PATH` environment variable.

This shouldn't depend on any external libraries.  Tested on CentOS 7 and Arch Linux (as of February 2024).

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
| -lt | Mean transition lifetime (in ps, use -ltns / -ltus / -lts / -lth for nanoseconds / microseconds / seconds / hours, respectively). |
| -hl | Transition half-life (in ps, use -hlns / -hlus / -hls / -hlh for nanoseconds / microseconds / seconds / hours, respectively) |
| -b | Reduced transition probability (for the L multipole) in units of e^2 fm^(2L) for electric multipoles or uN^2 fm^(2L-2) for magnetic multipoles. |

Optional parameters:

|**Parameter**|**Description**|
|:---:|:---:|
| -br | branching fraction of this transition (maximum 1, default 1) |
| -d | mixing ratio with the L+1 multipole |
| -ji | inital spin (integer or half-integer) |
| -jf | final spin (integer or half-integer) |
| -A | mass number of the nucleus |
| -Z | proton number of the nucleus |

Flags:

|**Flag**|**Description**|
|:---:|:---:|
| --barn | Use/calculate transition probability with spatial dimension in barns rather than fm (eg. B(E2) in e<sup>2</sup> b<sup>2</sup> rather than e<sup>2</sup> fm<sup>4</sup>). |
| --wu  |  Use/calculate transition probability in Weisskopf units (W.u.) rather than the default units specified above.  If used, requires the `-A` parameter. |
| --up | Use/calculate transition probability from final to initial state instead of vice versa.  If used, requires the `-ji` and `-jf` parameters. |
| --brrel | Specifies that the branching fraction provided with the `-br` option is actually an intensity relative to another transition. |
| --beta2 | Calculate the quadrupole deformation parameter, assuming a 2->0 (g.s.) transition.  Requires `-m E2 -ji 2 -jf 0`, and the `-A` and `-Z` parameters.  Assumes mean charge radius R = r_0*A^(1/3), with r_0 = 1.2 fm. |
| --quiet | Only show the result of the calculation. |
| --help | Print a list of parameters. |
