# JPerpleX

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://sc-dyer.github.io/JPerpleX.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://sc-dyer.github.io/JPerpleX.jl/dev/)
[![Build Status](https://github.com/sc-dyer/JPerpleX.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/sc-dyer/JPerpleX.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/sc-dyer/JPerpleX.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/sc-dyer/JPerpleX.jl)

This is a project still very much in development. The goal is to create an effective Perple_X wrapper in Julia. Currently, you can call meemum and do some pseudosection and isopleth plots but using this package requires local compilation of an unofficial fork of Perple_X that I have on my github: https://github.com/sc-dyer/Perple_X/tree/perplex-wrapper. To use this package right now you have to compile that specific branch (```make all``` or ```make perplexwrap```) and move the compiled binary ```perplexwrap.so``` to the local src directory in JPerplex.jl. I intend to submit a pull request to the official Perple_X project when this is a bit more mature so that the necessary functions can be registered to the Julia package manager and this local compilation step will be made unnecessary.
