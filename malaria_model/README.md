# malaria_model (draft source doc)

## To Build

`g++` is the default compiler. Several make targets are specified in the `Makefile`, but currently only the `release-cuda` and `test-cuda` targets can be compiled without modifications to the `Makefile`.


```
make
```
For more details please refer to the `Makefile`.

## To Run

```
./build/bin/malaria_model -c config.json -o outputs
```

### Mandatory Commandline Options

- `-c [FILE]`, user specified simulation configuration file. An example `config.json` is given below. See schema for more details.

```
{
    "village" :
    {
        "count" : 1000,
        "data_file" : 
        {
            "name" : "data/regions_data_demo.txt",
            "delimiter" : "t"
        }
    },
    "human" :
    {
        "male_percentage" : 0.48,
        "age_weight_vector_file": "data/age_weights_vector_demo.txt"
    }
}
```

- `-o [OUTPUT_DIRECTORY]`

### Optional Commandline Options

- `-i [OUTPUT_PREFIX]`, if given will be used as prefix for all outputs of this run, otherwise a timestamp will be used to identify outputs. 

## To Test
We use the [Catch2](https://github.com/catchorg/Catch2) framework.

1. Build tests:
```
make test
```

2. Run all tests:
```
malaria_model$ ./build/bin/test
```

3. List all tests:
```
malaria_model$ ./build/bin/test -l
```

4. Run a specific test:
```
malaria_model$ ./build/bin/test [blood_system:dominance_density_based]
```

## Dependencies
C++ Libraries and the version installed on the development machine:

- [Boost](http://www.boost.org/) Version 1.58 +

### .hpp Dependencies
These are header-only libraries that are included as part of this project. No installation is required.

- [Catch2](https://github.com/catchorg/Catch2), v2.1.2, a C++ test framework. `./test/third_party/`
- [RapidJSON](https://github.com/Tencent/rapidjson), v1.1.0, a fast JSON parser/generator for C++ with both SAX/DOM style API. `./include/third_party/`


## Next

- Age<->Mobility Relation
- Age by day i.e. distributed birthdays
- Stochasticity in Seasonality per village
- Probabilities such as 1/14 for incubation is too flat stochasticity