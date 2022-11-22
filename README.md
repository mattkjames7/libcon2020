# libcon2020

C++ implementation of the Connerney et al., 1981 and Connerney et al., 2020 Jovian magnetodisc model. This model provides the magnetic field due to a "washer-shaped" current disc near to Jupiter's magnetic equator. The model uses the analytical equations from Edwards et al., 2001 or Connerney et al., 1981; or the numerical integration of the Connerney et al., 1981 equations to provide the magnetic field vectors due to the azimuthal current. This code also implements the Connerney et al., 2020 radial current and the Leicester magnetosphere-ionosphere coupling (L-MIC) models which provide the azimuthal component of the mangetic field.

This is part of [libjupitermag](https://github.com/mattkjames7/libjupitermag.git), which is part of a greater effort to provide community code for the Jovian magnetosphere:

[Magnetospheres of the Outer Planets Group Community Code](https://lasp.colorado.edu/home/mop/missions/juno/community-code/)

## Building libcon2020

To build this library in Linux or Mac OS:

```bash
#clone this repo
git clone https://github.com/mattkjames7/libcon2020.git
cd libcon2020

#build
make 

#optionally install it system wide
sudo make install
```

In Windows:

```powershell
git clone https://github.com/mattkjames7/libcon2020.git
cd libcon2020

.\compile.bat
```

With a system wide installation, the compiler and linker will be able to locate the library and its header, otherwise absolute paths must be provided for linking and including. In Windows there is an experimental script ```install.bat``` which will copy the DLL and headers to folders within the `C:\TDC-GCC-64\` directory. This is experimental, instead it might be better to copy the headers and the DLL to the root directory of the executable linked to it.

Uninstallation can be acheived in Linux and Mac using ```sudo make uninstall```.

## Usage

### Linking to and Including libcon2020

If a system-wide installation is successful then the library may be linked to simply by including the ```-lcon2020``` flag while compiling/linking. Otherwise the path to the library must also be included, e.g. ```-L /path/to/lib/directory -lcon2020```. In Windows, the DLL should be placed in the root directory of the executable linked to it.

This library includes two headers: `con2020.h` - a full header which provides access to everything within the library for use with C++ (including classes and C++ specific objects); and `con2020c.h` - a partial header which can be included in C, exposing C-compatible wrapper functions. The wrapper functions in the partial header file would also provide the easiest ways to link other languages to the library such as Python, IDL and Fortran.

If the library was installed system-wide, then the headers may be included using ```#include <con2020.h>``` and ```#include <con2020c.h``` for C++ and C, respectively. Otherwise, a relative or absolute path to the headers must be used, e.g. ```#include "path/to/con2020.h"```.

### C++ usage

This section briefly describes some C++ specific examples for using the `libcon2020` library, while the following section is also somewhat applicable.

### Other Languages
