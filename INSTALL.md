# Installation

We explain two different ways to install *cubic_interpolation*: With or without the package manager [*conan*](https://conan.io/).

### Approach 1: Installation with conan

For this installation approach, all dependencies will be fetched by the package manager conan.
This way, you don't have to install them by yourself.

If conan is not yet installed, you can do so, for example by using pip:

```
$ pip install conan
```

Clone the repository and create a build directory

```
$ git clone https://github.com/MaxSac/cubic_interpolation
$ cd cubic_interpolation && mkdir build && cd build      
```
Use conan to fetch all dependencies:

```
$ conan install ..
```

Next, build and install *cubic_interpolation*:

```
$ conan build ..
# cmake --install .
```

### Approach 2: Installation without conan

If you do not want to use conan, you need to provide these dependencies:

* [boost](https://www.boost.org/) (>1.75.0)
* [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) (>3.3.9)

If you have these dependencies installed, clone the repository and create a build directory:

```
$ git clone https://github.com/MaxSac/cubic_interpolation
$ cd cubic_interpolation && mkdir build && cd build      
```

Use *CMake* and *make* to build and install *cubic_interpolation*:

```
$ cmake .. -DCMAKE_BUILD_TYPE=Release	# or other CMake options
$ make -j
# make install
```