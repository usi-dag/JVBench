## JVBench Benchmark Suite

JVBench is a collection of diverse benchmarks
that exercise the [Java Vector API](https://openjdk.org/jeps/426) on typical SIMD workloads.

This benchmark suite includes a version of the workloads of the 
[RIVEC](https://github.com/RALC88/riscv-vectorized-benchmark-suite) suite after
recasting them in Java and expressing vector operations with the Java Vector API.


### List of benchmarks

The following table is a complete list of benchmarks included in JVBench in alphabetical order.

| Application Name | Application Domain         | Algorithmic Model    | Taken From |
|------------------|----------------------------|----------------------|------------|
| Axpy             | High Performance Computing | BLAS                 | RIVEC      |
| Blackscholes     | Financial Analysis         | Dense Linear Algebra | PARSEC     |
| Canneal          | Engineering                | Unstructured Grids   | PARSEC     |
| Jacobi-2D        | Engineering                | Dense Linear Algebra | PolyBench  |
| LavaMD           | Molecular Dynamics         | N-Body               | Rodinia    |
| Particlefilter   | Medical Imaging            | Structured Grids     | Rodinia    |
| Pathfinder       | Grid Traversal             | Dynamic Programming  | Rodinia    |
| Somier           | Physics Simulation         | Dense Linear Algebra | RIVEC      |
| Streamcluster    | Data Mining                | Dense Linear Algebra | PARSEC     |
| Swaptions        | Financial Analysis         | MapReduce Regular    | PARSEC     |


### Prerequisites
- [OpenJDK HotSpot JDK 16](https://github.com/openjdk/jdk16/) or above


### Obtaining the suite
To run the suite, you can download the suite JAR [here](http://195.176.181.79/cgo23/JVBench-1.0.jar).
Alternatively, you can build it yourself. To do, you will need [Maven](https://maven.apache.org/) version 3.8 or above. 
To build an executable JAR,  use the following command:
```shell
$ mvn clean package
```


### Running the suite
To run a JVBench benchmark, execute the following java command:
```shell
$ java --add-modules jdk.incubator.vector -jar JVBench-1.0.jar "<benchmarks>"
```
where <benchmarks> is the benchmark name that you wish to run. Please append the word "Benchmark" to the benchmark name. For example, to run benchmark Axpy, specify `AxpyBenchmark` as the benchmark.

By default, the suite executes each benchmark operation for a specific number of times.
For thorough experimental evaluation, the benchmarks should be repeated for a large number of times or executed for a long time. 
The number of repetitions and the execution time can be set for all benchmarks using the -f, -wi, and -i options.
For a complete description refers to [JHM tutorial](https://github.com/guozheng/jmh-tutorial/blob/master/README.md)

Moreover, it is possible to override the default benchmark input with the following command.

```shell
$ java --add-modules jdk.incubator.vector <inputs> -jar JVBench-1.0.jar "<benchmarks>"
```

For example, to specify an input size of 70k elements to Axpy, you can run the following command:

```shell
$ java --add-modules -Dsize=70000 jdk.incubator.vector -jar JVBench-1.0.jar "AxpyBenchmark"
```

The default input (expressed as different system properties) for each benchmark is listed in the table below

| Application Name | Default Input                                                                                                                         |
|------------------|---------------------------------------------------------------------------------------------------------------------------------------|
| Axpy             | -Dsize=70000                                                                                                                          |
| Blackscholes     | -Dinput=blackscholes/input/in_512K.input                                                                         |
| Canneal          | -Dnswaps=10000 -DTEMP=2000 -Dnetlist=canneal/input/2500000.nets -Dnsteps=300                                     |
| Jacobi-2D        | -Dsize=10000 -Dtsteps=14                                                                                                              |
| LavaMD           | -Dinput=lavaMD/input/lavaMD_127776.input                                                                           |
| Particlefilter   | -Dx=128 -Dy=128 -Dz=24 -Dnp=32768                                                                                                     |
| Pathfinder       | -Dinput=pathfinder/input/pathfinder_5000_5000.input                                                                |
| Somier           | -Dsteps=10 -Dn=128                                                                                                                    |
| Streamcluster    | -Dk1=3 -Dk2=10 -Ddim=128 -Dchunksize=128 -Dclustersize=10 -Dinput=streamcluster/input/streamcluster_128_128.input |
| Swaptions        | -Dms=64 -Dns=16384                                                                                                                    |


### 🤝 Contributing
You are welcome to contribute to this open-source benchmark suite. If you want to add or modify
the code, feel free to fork the repo, create issues and finally a pull request to review your code!

