[![DOI](https://zenodo.org/badge/33890374.svg)](https://zenodo.org/badge/latestdoi/33890374)

# trackcpp

Particle tracking code

## Code check

In order to use clang-tidy, a compile_commands.json file must be present at the root of the repository.
In order to generate one and still use make as the build system:

```command
# Compile using 'bear <build command>'
bear make -j$(nproc)
```
