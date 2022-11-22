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

compile.bat
```

With a system wide installation, the compiler and linker will be able to locate the library and its header, otherwise absolute paths must be provided for linking and including.
