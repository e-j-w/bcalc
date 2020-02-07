# bcalc

Reduced transition probability calculator

Maintainer: Jonathan Williams

## Installation

Use `make` to compile. To run the program from anywhere, move the resulting `bcalc` executable to any directory under your `$PATH` environment variable.

This shouldn't depend on any external libraries.  Tested on CentOS 7.

## Usage

The calculator is run from the command line:

```
bcalc -E VALUE -M VALUE -Lt VALUE -Hl VALUE -B VALUE
```
### Parameters

Both of the following are needed:

|**Parameter**|**Description**|
|:---:|:---:|
| -E | transition energy in keV |
| -M |  multipole (eg. E1, M1, E2, etc.) |

One of the following is needed:

|**Parameter**|**Description**|
|:---:|:---:|
| -Lt | Mean transition lifetime (in ps) |
| -Hl | Transition half-life (in ps) |
| -B | Reduced transition probability |

Optional parameters:

|**Parameter**|**Description**|
|:---:|:---:|
| -up | Use/calculate transition probability from final to initial state instead of vice versa.  If used, requires the `-ji` and `-jf` parameters. |
| -verbose | Print extra info |
| -ji | inital spin |
| -jf | final spin |

Mixed transitions are not yet implemented.