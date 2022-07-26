# C++ Project Template

A simple template repository to start a new C++ project using CMake.

Click on the green [`Use this template`](https://github.com/ssciwr/cpp-project-template/generate) button to get started.

If you are looking for more advanced features (such as Python bindings or integration with sites like ReadTheDocs, codecov, sonarcloud or PyPI)
take a look at our [C++ Project Cookiecutter](https://github.com/ssciwr/cookiecutter-cpp-project)

## Contents

This example project contains the `adder` library,
an application `adder_app` which uses this library,
and a test-suite which tests the library.

Any pull-requests or commits to the repository trigger GitHub Actions,
which will compile the code and run the tests.

Project structure:

- [src](src)
  - the `adder` library source code
  - this is where the meat of the project is: the implementation
- [include/adder](include/adder)
  - the `adder` library headers
  - the public interface of the library
- [app](app)
  - the application which uses the `adder` library
- [tests](tests)
  - the test code
  - each `x.cpp` file has a corresponding `x_t.cpp` file here with tests
- [ext](ext)
  - external libraries, e.g. Catch2 testing framework
- [.github/workflows/ci.yml](.github/workflows/ci.yml)
  - the GitHub Actions configuration

## Compiling

To compile the project and run the tests:

```
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make
make test
```

## Documentation

If you have Doxygen installed you can also build the documentation by enabling the `BUILD_DOCS` CMake option, and then running `make doxygen`:

```
cmake .. -DCMAKE_BUILD_TYPE=Release -DBUILD_DOCS=ON
make doxygen
```

This will generate the documentation in the `html` folder.