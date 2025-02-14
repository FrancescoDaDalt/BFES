# What is this Repository

This repository contains all code related to the publication of the BFES algorithm.
The primarly language is C++ since it yields the highest performance and flexibility under reasonable development time. It furthermore ensures an even playfield for the various benchmarked algorithms since they vary significantly in complexity.
A secondary implementation of the BFES insertion routine is provided in Java for interoperation with Apache Flink.

Regarding licenses, see each file. Broadly speaking, all code written by the main author is under MIT-license with exception of other people's algorithms (see references) to which some other license may apply.

# Prerequisites

Download this repository and navigate to the directory containing this file.

# Running C++ benchmarks

## 1. Dependencies

Install Cmake, GSL, Eigen3, OpenMP.

Install MOSEK if possible. MOSEK is required for running BFES with "u>1" and SeqSketch.
MOSEK requires a license but licenses are free for academic purposes.
If MOSEK is not installed, the code will still compile but omit SeqSketch and do a compile-time check that u==1.

## 2. Mac specific

Change the line ```set(LIBOMP_DIR "/opt/homebrew/Cellar/libomp/19.1.5")``` in the CMakeLists.txt file to point to your OpenMP installation.

## 3. MOSEK specific

Change the line ```set(MOSEK_DIR "/opt/mosek/10.0/tools/platform/arch")``` to point to your MOSEK installation directory.

## 4. Configure

Modify the ```main_benchmarks.cpp``` file to your desire. By default it synthesizes a datastream with a Gaussian distribution and 1M unique keys. Modify the ```U``` in ```using hypergrid_type = hypergrid<U>;``` to adjust the "u" parameter from BFES. If MOSEK is not available, u has to be 1.

## 5. Compile

```mkdir -p build```

```cd build```

```cmake ..```

```make```

## 6. MOSEK + Mac specific

Inside the build directory, run ```export DYLD_LIBRARY_PATH=<PATH_TO_MOSEK>/bin``` where ```<PATH_TO_MOSEK>``` is the path from point 3. 

## Run the executable

```./BFES```

# Running in Docker

## 1. Requirements

Install Docker and start the docker deamon.

## 2. Build the Docker image, start container, and run benchmarks

```docker build -t bfes_image .```

```docker run -it bfes_image```

This gives you a shell in ```/app/build``` where the ```C++``` executable is compiled.
You can run it by calling ```./BFES```.
Docker will fail if "u>1" because the Dockerfile does not install MOSEK. In that case, adjust the parameter and try again.

# Apache Flink

We have added a ```ProcessFunction``` implementation of the BFES insertion routine. 
It processes key-value tuples and populates the BFES sketch datastructure.
The output of the stream is a stream of sketch datastructure snapshots which contain 10 ms (can be adjusted based on necessities) worth of stream information.
To compute the estimates, the sketch snapshots must be passed to the ```C++``` implmentation.

## Run Flink in Docker

```mvn package -f /app/build/../ApacheFlink_x_Java/pom.xml```

```/opt/flink/bin/start-cluster.sh```

```/opt/flink/bin/flink run /app/ApacheFlink_x_Java/target/FlinkBFES-job-1.0-SNAPSHOT-jar-with-dependencies.jar```

# References
The following references list all non-standard algorithms / libraries that were used either directly or indirectly for benchmarking:

1. A. Appleby, "MurmurHash3," SMHasher. [Online]. Available: https://github.com/aappleby/smhasher/blob/master/src/MurmurHash3.cpp. [Accessed: 14-Feb-2025].

2. G. Guennebaud, B. Jacob, et al, "Eigen: A C++ template library for linear algebra." [Online]. Available: https://eigen.tuxfamily.org. [Accessed: 14-Feb-2025].

3. M. Galassi et al, "GNU Scientific Library (GSL)," Free Software Foundation. [Online]. Available: https://www.gnu.org/software/gsl/. [Accessed: 14-Feb-2025].

4. MOSEK ApS, "The MOSEK optimization toolbox." [Online]. Available: https://www.mosek.com. [Accessed: 14-Feb-2025].

5. Apache Software Foundation, "Apache Flink: Stream processing framework." [Online]. Available: https://flink.apache.org. [Accessed: 14-Feb-2025].

6. Apache Software Foundation, "Apache DataSketches: Scalable approximation algorithms." [Online]. Available: https://datasketches.apache.org. [Accessed: 14-Feb-2025].

7. D. Cai, M. Mitzenmacher, and R. P. Adams, “A Bayesian nonparametric view on count-min sketch,” in Advances in Neural Information Processing Systems, vol. 31, S. Bengio, H. Wallach, H. Larochelle, K. Grauman, N. Cesa-Bianchi, and R. Garnett, Eds. Curran Associates, Inc., 2018. [Online]. Available: https://proceedings.neurips.cc/paper_files/paper/2018/file/0b9e57c46de934cee33b0e8d1839bfc2-Paper.pdf.

8. T. Akiba, “HLL-HIP,” 2014. [Online]. Available: https://github.com/iwiwi/hyperloglog-hip. Accessed: Jun. 10, 2024.

9. M. Charikar, K. Chen, and M. Farach-Colton, “Finding frequent items in data streams,” in Automata, Languages and Programming, Berlin, Heidelberg: Springer Berlin Heidelberg, 2002, pp. 693–703.

10. G. Cormode and S. Muthukrishnan, “An improved data stream summary: The count-min sketch and its applications,” J. Algorithms, vol. 55, no. 1, pp. 58–75, 2005. doi: 10.1016/j.jalgor.2003.12.001.

11. P. Flajolet, E. Fusy, O. Gandouet, and F. Meunier, “HyperLogLog: The analysis of a near-optimal cardinality estimation algorithm,” Discrete Math. Theor. Comput. Sci., vol. AH, Mar. 2012, doi: 10.46298/dmtcs.3545.

12. R. Gallager, “Low-density parity-check codes,” IRE Trans. Inf. Theory, vol. 8, no. 1, pp. 21–28, 1962. doi: 10.1109/TIT.1962.1057683.

13. S. Geman and D. Geman, “Stochastic relaxation, Gibbs distributions, and the Bayesian restoration of images,” IEEE Trans. Pattern Anal. Mach. Intell., vol. PAMI-6, no. 6, pp. 721–741, 1984. doi: 10.1109/TPAMI.1984.4767596.

14. F. D. Da Dalt, S. Scherrer, and A. Perrig, “Bayesian sketches for volume estimation in data streams,” Proc. VLDB Endow., vol. 16, no. 4, pp. 657–669, Dec. 2022. doi: 10.14778/3574245.3574252.

15. S. Sheng, Q. Huang, S. Wang, and Y. Bao, “PR-Sketch: Monitoring per-key aggregation of streaming data with nearly full accuracy,” Proc. VLDB Endow., vol. 14, no. 10, pp. 1783–1796, Jun. 2021. doi: 10.14778/3467861.3467868.

16. Q. Huang, S. Sheng, X. Chen, Y. Bao, R. Zhang, Y. Xu, and G. Zhang, “Toward nearly-zero-error sketching via compressive sensing,” in Proc. 18th USENIX Symp. Networked Syst. Des. Implement. (NSDI '21), 2021, pp. 1027–1044. [Online]. Available: https://www.usenix.org/conference/nsdi21/presentation/huang.

17. E. Dolera, S. Favaro, and S. Peluchetti, “A Bayesian nonparametric approach to count-min sketch under power-law data streams,” arXiv preprint, arXiv:2102.03743, 2021. [Online]. Available: https://arxiv.org/abs/2102.03743.