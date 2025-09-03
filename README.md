# ECG2O
**Equality-Constrained g2o – Extending g2o with Equality Constraint Support**

ECG2O extends the open-source **g2o (General Graph Optimization)** framework by introducing **equality constraints**. 
This allows for the addition of **equality edge factors**, enabling constrained optimization within factor graphs. 
The library also includes an **Augmented Lagrangian implementation** as a benchmark for constrained optimization techniques.

### Reference

> *ecg2o: A Seamless Extension of g2o for Equality-Constrained Factor Graph Optimization*  

---

## 📥 Clone the Repository

```bash
git clone https://github.com/snt-arg/ecg2o.git
cd ecg2o
```

## 🚀 Quickstart with Docker

### 1. Build the Docker image

```bash
docker build -t ecg2o -f docker/Dockerfile .
```

### 2. Run the container

```bash
docker run -it --rm ecg2o
```

### 3. Inside the container

```bash
cd /workspace
cmake .. 
make
```


### 4. Run the examples
```bash
./example_gn
```
or
```bash
./example_al
```
or
```bash
./example_soft
```

or
```bash
./oc_test
```


---

## 🔧 Manual Installation Without Docker

To set up the environment without Docker, follow these steps:

1. **Install Dependencies (Ubuntu 22.04)**

   ```bash
   sudo apt update
   sudo apt install -y \
       git cmake build-essential libeigen3-dev libspdlog-dev \
       libsuitesparse-dev qtdeclarative5-dev qt5-qmake \
       libqglviewer-dev-qt5 libmetis-dev
   ```

2. **Clone and Build g2o with Patch**

   ```bash
   cd /opt
   git clone https://github.com/RainerKuemmerle/g2o
   cd g2o
   git checkout master
   ```

3. **Apply Required Code Patch**

   Edit the file `g2o/core/sparse_optimizer.h`:

   Find the line:
   ```cpp
     void update(const double* update);
   ```

   Replace it with:
   ```cpp
     virtual void update(const double* update);
   ```

   This change enables overriding the `update` method required by the BIPM extension.

4. **Build and Install g2o**

   ```bash
   mkdir build
   cd build
   cmake -DCMAKE_BUILD_TYPE=Release -DG2O_USE_OPENMP=ON -DCMAKE_CXX_FLAGS="-DEIGEN_USE_THREADS -DEIGEN_USE_OPENMP" ..
   make -j$(nproc)
   sudo make install
   ```

5. **Build Your Application**

   Assuming you're in your app's build directory:

   ```bash
   cmake ..
   make
    ```




6. **Usage

Once the library is successfully built, you can execute the provided functions for equality-constrained optimization.

```bash
./example_gn
```
or
```bash
./example_al
```
or
```bash
./example_soft
```

or
```bash
./oc_test
```

Ensure that all necessary input files or configurations are correctly set up before execution.

---



## License

This project is open-source and distributed under GPL-3.0
.
