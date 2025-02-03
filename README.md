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
It processes key-value tuples and poplates the BFES sketch datastructure.
The output of the stream is a stream of sketch datastructure snapshots which contain 10 ms (can be adjusted based on necessities) worth of stream information.
To compute the estimates, the sketch snapshots must be passed to the ```C++``` implmentation.

## Run Flink in Docker

```mvn package -f /app/build/../ApacheFlink_x_Java/pom.xml```

```/opt/flink/bin/start-cluster.sh```

```/opt/flink/bin/flink run /app/ApacheFlink_x_Java/target/FlinkBFES-job-1.0-SNAPSHOT-jar-with-dependencies.jar```

